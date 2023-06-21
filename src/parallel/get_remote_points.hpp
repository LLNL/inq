/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__GET_REMOTE_POINTS
#define INQ__PARALLEL__GET_REMOTE_POINTS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/run.hpp>
#include <gpu/array.hpp>
#include <parallel/communicator.hpp>
#include <parallel/array_iterator_2d.hpp>
#include <parallel/partition.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <cstdlib> //drand48

namespace inq {

namespace basis {
template<class BasisType, typename ElementType> class field;
template<class BasisType, typename ElementType, class PartitionType> class field_set;
}

namespace parallel {

struct remote_points_table {

	struct point_position {
		long point;
		long position;
	};
	
	std::vector<std::vector<point_position>> points_needed;
	int total_needed;
	int total_requested;
	gpu::array<int, 1> list_points_requested;
	gpu::array<int, 1> list_sizes_requested;
	gpu::array<int, 1> list_sizes_needed;	
	gpu::array<int, 1> list_needed;
	gpu::array<int, 1> list_displs_needed;
	gpu::array<int, 1> list_displs_requested;	
	
	template <class BasisType, class ArrayType>
	remote_points_table(BasisType & basis, ArrayType const & point_list){

		CALI_CXX_MARK_FUNCTION;

		auto const num_proc = basis.comm().size();

		// create a list per processor of the points we need
		points_needed.resize(num_proc);
		
		for(long ilist = 0; ilist < long(point_list.size()); ilist++){
			auto point = point_list[ilist];
			assert(point >= 0);
			assert(point < basis.size());

			auto src_proc = basis.part().location(point);
			assert(src_proc >= 0 and src_proc < num_proc); 
			points_needed[src_proc].push_back(point_position{point, ilist});
		}

		// find out how many points each processor requests from us
		list_sizes_needed.reextent(num_proc);
		list_sizes_requested.reextent(num_proc);	
	
		total_needed = 0;
		for(int iproc = 0; iproc < num_proc; iproc++){
			list_sizes_needed[iproc] = points_needed[iproc].size();
			total_needed += list_sizes_needed[iproc];
		}

		assert(total_needed == long(point_list.size()));

		MPI_Alltoall(raw_pointer_cast(list_sizes_needed.data_elements()), 1, MPI_INT, raw_pointer_cast(list_sizes_requested.data_elements()), 1, MPI_INT, basis.comm().get());

		// get the list of points each processor requests from us and send a list of the points we need
		list_needed.reextent(total_needed);
		list_displs_needed.reextent(num_proc);
		list_displs_requested.reextent(num_proc);
	
		total_requested = 0;

		basis.comm().barrier();
	
		for(int iproc = 0; iproc < num_proc; iproc++){

			total_requested += list_sizes_requested[iproc];
		
			if(iproc > 0){
				list_displs_needed[iproc] = list_displs_needed[iproc - 1] + list_sizes_needed[iproc - 1];
				list_displs_requested[iproc] = list_displs_requested[iproc - 1] + list_sizes_requested[iproc - 1];			
			} else {
				list_displs_needed[iproc] = 0;
				list_displs_requested[iproc] = 0;
			}

		}

		basis.comm().barrier();
	
		for(int iproc = 0; iproc < num_proc; iproc++){

			basis.comm().barrier();
		
			for(long ip = 0; ip < list_sizes_needed[iproc]; ip++) {
				list_needed[list_displs_needed[iproc] + ip] = points_needed[iproc][ip].point;
			}
		}

		list_points_requested.reextent(total_requested);
	
		MPI_Alltoallv(raw_pointer_cast(list_needed.data_elements()), raw_pointer_cast(list_sizes_needed.data_elements()), raw_pointer_cast(list_displs_needed.data_elements()), MPI_INT,
									raw_pointer_cast(list_points_requested.data_elements()), raw_pointer_cast(list_sizes_requested.data_elements()), raw_pointer_cast(list_displs_requested.data_elements()), MPI_INT, basis.comm().get());
		
	}

};

template <class BasisType, class ElementType, class ArrayType>
gpu::array<ElementType, 1> get_remote_points(basis::field<BasisType, ElementType> const & source, ArrayType const & point_list){

	CALI_CXX_MARK_FUNCTION;

	auto const num_proc = source.basis().comm().size();
 
	if(num_proc == 1){
		gpu::array<ElementType, 1> remote_points(point_list.size());
		
		gpu::run(point_list.size(),
						 [rem = begin(remote_points), sou = begin(source.linear()), poi = begin(point_list)] GPU_LAMBDA (auto ip){
							 rem[ip] = sou[poi[ip]];
						 });
		
		return remote_points;
	}
	
	remote_points_table rp(source.basis(), point_list);
	
	// send the value of the request points
	gpu::array<ElementType, 1> value_points_requested(rp.total_requested);
	gpu::array<ElementType, 1> value_points_needed(rp.total_needed);

	gpu::run(rp.total_requested,
					 [par =  source.basis().part(), lreq = begin(rp.list_points_requested), vreq = begin(value_points_requested), sou = begin(source.linear())] GPU_LAMBDA (auto ip){
						 auto iplocal = par.global_to_local(parallel::global_index(lreq[ip]));
						 vreq[ip] = sou[iplocal];
					 });
	
	auto mpi_type = boost::mpi3::detail::basic_datatype<ElementType>();
	
	MPI_Alltoallv(raw_pointer_cast(value_points_requested.data_elements()), raw_pointer_cast(rp.list_sizes_requested.data_elements()), raw_pointer_cast(rp.list_displs_requested.data_elements()), mpi_type,
								raw_pointer_cast(value_points_needed.data_elements()), raw_pointer_cast(rp.list_sizes_needed.data_elements()), raw_pointer_cast(rp.list_displs_needed.data_elements()), mpi_type, source.basis().comm().get());

	// Finally copy the values to the return array in the proper order
	gpu::array<ElementType, 1> remote_points(point_list.size());

	long ip = 0;
	for(int iproc = 0; iproc < num_proc; iproc++){
		for(long jp = 0; jp < long(rp.points_needed[iproc].size()); jp++) {
			remote_points[rp.points_needed[iproc][jp].position] = value_points_needed[ip];
			ip++;
		}
	}

	assert(ip == long(point_list.size()));
	
  return remote_points;
}

template <typename FieldSetType, typename ArrayType>
gpu::array<typename FieldSetType::element_type, 2> get_remote_points(FieldSetType const & source, ArrayType const & point_list){

	CALI_CXX_MARK_FUNCTION;

	auto const nset = source.local_set_size();
	auto const num_proc = source.basis().comm().size();
	
	if(num_proc == 1){
		gpu::array<typename FieldSetType::element_type, 2> remote_points({point_list.size(), nset});
		
		gpu::run(nset, point_list.size(),
						 [rem = begin(remote_points), sou = begin(source.matrix()), poi = begin(point_list)] GPU_LAMBDA (auto ist, auto ip){
							 rem[ip][ist] = sou[poi[ip]][ist];
						 });
		
		return remote_points;
	}
	
	remote_points_table rp(source.basis(), point_list);
	
	// send the value of the request points
	gpu::array<typename FieldSetType::element_type, 2> value_points_requested({rp.total_requested, nset});
	gpu::array<typename FieldSetType::element_type, 2> value_points_needed({rp.total_needed, nset});

	gpu::run(nset, rp.total_requested,
					 [par =  source.basis().part(), lreq = begin(rp.list_points_requested), vreq = begin(value_points_requested), sou = begin(source.matrix())] GPU_LAMBDA (auto iset, auto ip){
						 auto iplocal = par.global_to_local(parallel::global_index(lreq[ip]));
						 vreq[ip][iset] = sou[iplocal][iset];
					 });
	
	MPI_Datatype mpi_type;
	MPI_Type_contiguous(nset, boost::mpi3::detail::basic_datatype<typename FieldSetType::element_type>(), &mpi_type);
	MPI_Type_commit(&mpi_type);
	
	MPI_Alltoallv(raw_pointer_cast(value_points_requested.data_elements()), raw_pointer_cast(rp.list_sizes_requested.data_elements()), raw_pointer_cast(rp.list_displs_requested.data_elements()), mpi_type,
								raw_pointer_cast(value_points_needed.data_elements()), raw_pointer_cast(rp.list_sizes_needed.data_elements()), raw_pointer_cast(rp.list_displs_needed.data_elements()), mpi_type, source.basis().comm().get());

	MPI_Type_free(&mpi_type);
	
	// Finally copy the values to the return array in the proper order
	gpu::array<typename FieldSetType::element_type, 2> remote_points({point_list.size(), nset});

	long ip = 0;
	for(int iproc = 0; iproc < num_proc; iproc++){
		for(long jp = 0; jp < long(rp.points_needed[iproc].size()); jp++) {
				for(long iset = 0; iset < nset; iset++) remote_points[rp.points_needed[iproc][jp].position][iset] = value_points_needed[ip][iset];
			ip++;
		}
	}

	assert(ip == long(point_list.size()));
	
  return remote_points;
}

template <typename FieldSetType, typename ArrayType>
gpu::array<typename FieldSetType::element_type, 2> get_remote_points(FieldSetType const & source, ArrayType const & point_list, ArrayType const & state_list){

	CALI_CXX_MARK_FUNCTION;

	gpu::array<typename FieldSetType::element_type, 2> remote_points({point_list.size(), state_list.size()});
	
	if(source.full_comm().size() == 1) {
		gpu::run(state_list.size(), point_list.size(),
						 [rem = begin(remote_points), sou = begin(source.matrix()), poi = begin(point_list), sta = begin(state_list)] GPU_LAMBDA (auto ist, auto ip){
							 rem[ip][ist] = sou[poi[ip]][sta[ist]];
						 });

	} else {

		for(parallel::array_iterator_2d pai(source.basis().part(), source.set_part(), source.full_comm(), source.matrix()); pai != pai.end(); ++pai){

			for(long ip = 0; ip < point_list.size(); ip++){
				if(not source.basis().part().contains(point_list[ip], pai.xpart())) continue;

				for(long ist = 0; ist < state_list.size(); ist++){
					if(not source.set_part().contains(state_list[ist], pai.ypart())) continue;
					
					remote_points[ip][ist] = (*pai)[point_list[ip] - source.basis().part().start(pai.xpart())][state_list[ist] - source.set_part().start(pai.ypart())];
				}
      }
			
		}
	}
		
	return remote_points;
}

}
}
#endif

#ifdef INQ_PARALLEL_GET_REMOTE_POINTS_UNIT_TEST
#undef INQ_PARALLEL_GET_REMOTE_POINTS_UNIT_TEST

#include <math/vector3.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>

#include <catch2/catch_all.hpp>
#include <parallel/communicator.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;
	
  parallel::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	auto set_comm = basis::set_subcomm(cart_comm);
	auto basis_comm = basis::basis_subcomm(cart_comm);	

	systems::box box = systems::box::orthorhombic(13.3_b, 6.55_b, 8.02_b);
  basis::real_space rs(box, /*spacing =*/ 0.49672941, cart_comm);

	int const nvec = 19;
	
  basis::field<basis::real_space, complex> test_field(rs);
  basis::field_set<basis::real_space, complex> test_field_set(rs, 19);

  for(long ip = 0; ip < rs.local_size(); ip++){
    auto ipg = rs.part().local_to_global(ip);
    test_field.linear()[ip] = complex(ipg.value(), 0.1*ipg.value());
    for(int ivec = 0; ivec < nvec; ivec++) test_field_set.matrix()[ip][ivec] = (ivec + 1.0)*complex(ipg.value(), 0.1*ipg.value());		
  }

	srand48(500 + 34895783*cart_comm.rank());
	
	long const npoints = 1 + drand48()*(rs.size() - 1);
	assert(npoints > 0);
	assert(npoints <= rs.size());
	
	gpu::array<long, 1> list(npoints);
	
	for(long ip = 0; ip < npoints; ip++){
		list[ip] = drand48()*(rs.size() - 1);
		assert(list[ip] < rs.size());
	}
	
	SECTION("field"){
		auto remote_points = parallel::get_remote_points(test_field, list);
		
		for(long ip = 0; ip < npoints; ip++){
			CHECK(double(list[ip]) == Approx(real(remote_points[ip])));
			CHECK(0.1*double(list[ip]) == Approx(imag(remote_points[ip])));
		}
	}
	
	SECTION("field_set"){
		auto remote_points_set = parallel::get_remote_points(test_field_set, list);	

		for(long ip = 0; ip < npoints; ip++){
			for(int ivec = 0; ivec < nvec; ivec++){
				CHECK((ivec + 1.0)*double(list[ip]) == Approx(real(remote_points_set[ip][ivec])));
				CHECK((ivec + 1.0)*0.1*double(list[ip]) == Approx(imag(remote_points_set[ip][ivec])));
			}
		}
	}
	
	long const nst = 1 + drand48()*(nvec - 1);
	assert(nst > 0);
	assert(nst <= nvec);
		
	gpu::array<long, 1> stlist(nst);
	
	for(long ist = 0; ist < nst; ist++){
		stlist[ist] = drand48()*(nvec - 1);
		assert(stlist[ist] < nvec);
	}
	
	SECTION("field_set state indices"){
		
		auto remote_points_set = parallel::get_remote_points(test_field_set, list, stlist);	
		
		assert(remote_points_set.size() == npoints);
		assert(remote_points_set.transposed().size() == nst);

		for(long ip = 0; ip < npoints; ip++){
			for(long ist = 0; ist < nst; ist++){
				CHECK((stlist[ist] + 1.0)*double(list[ip]) == Approx(real(remote_points_set[ip][ist])));
				CHECK((stlist[ist] + 1.0)*0.1*double(list[ip]) == Approx(imag(remote_points_set[ip][ist])));
			}
		}
	}
	
}
#endif
