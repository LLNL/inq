/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__TRANSFER
#define OPERATIONS__TRANSFER

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <basis/fourier_space.hpp>
#include <basis/field_set.hpp>
#include <gpu/run.hpp>
#include <parallel/get_remote_points.hpp>
#include <operations/space.hpp>

#include <multi/adaptors/fftw.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/cufft.hpp>
#endif

#include <cassert>

namespace inq {
namespace operations {
namespace transfer {

template <class FieldType>
FieldType enlarge(FieldType const & source, typename FieldType::basis_type const & new_basis, double const factor = 1.0) {

	CALI_CXX_MARK_FUNCTION;
	
	FieldType destination(new_basis);
	destination.fill(0.0);
	
	if(not source.basis().part().parallel()){

		gpu::run(source.basis().sizes()[2], source.basis().sizes()[1], source.basis().sizes()[0],
						 [sbas = source.basis().point_op(), dbas = destination.basis().point_op(), des = begin(destination.cubic()), sou = begin(source.cubic()), factor] GPU_LAMBDA (auto iz, auto iy, auto ix){
							 
							 auto ii = sbas.to_symmetric_range(ix, iy, iz);
							 auto idest = dbas.from_symmetric_range(ii);
							 
							 des[idest[0]][idest[1]][idest[2]] = factor*sou[ix][iy][iz];
						 });
		
	} else {

		std::vector<long> point_list;

		for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < source.basis().sizes()[2]; iz++){
					
					auto ii = source.basis().to_symmetric_range(ix, iy, iz);
					auto idest = destination.basis().from_symmetric_range(ii);

					if(not destination.basis().local_contains(idest)) continue;

					point_list.push_back(source.basis().linear_index(ix, iy, iz));
					
				}
			}
		}

		gpu::array<long, 1> list(point_list.begin(), point_list.end());
		
		auto points = parallel::get_remote_points(source, list);
		
		long ip = 0;
		for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < source.basis().sizes()[2]; iz++){

					auto ii = source.basis().to_symmetric_range(ix, iy, iz);
					auto idest = destination.basis().from_symmetric_range(ii);

					if(not destination.basis().local_contains(idest)) continue;

					auto il0 = destination.basis().cubic_part(0).global_to_local(parallel::global_index(idest[0]));
					auto il1 = destination.basis().cubic_part(1).global_to_local(parallel::global_index(idest[1]));
					auto il2 = destination.basis().cubic_part(2).global_to_local(parallel::global_index(idest[2]));

					destination.cubic()[il0][il1][il2] = factor*points[ip];
					ip++;
					
				}
			}
		}

	}

	return destination;

}

//////////////////////////////////////////////////////////

template <class Type, class BasisType>
basis::field_set<BasisType, Type> enlarge(basis::field_set<BasisType, Type> const & source, BasisType const & new_basis, double const factor = 1.0) {
	CALI_CXX_MARK_FUNCTION;
	
	basis::field_set<BasisType, Type> destination(new_basis, source.set_size(), source.full_comm());
	destination.fill(0.0);

	if(not source.basis().part().parallel()){
		
		gpu::run(source.basis().sizes()[2], source.basis().sizes()[1], source.basis().sizes()[0],
						 [nst = source.set_part().local_size(), sbas = source.basis().point_op(), dbas = destination.basis().point_op(), des = begin(destination.hypercubic()), sou = begin(source.hypercubic()), factor]
						 GPU_LAMBDA (auto iz, auto iy, auto ix){
							 
							 auto ii = sbas.to_symmetric_range(ix, iy, iz);
							 auto idest = dbas.from_symmetric_range(ii);
							 
							 for(int ist = 0; ist < nst; ist++) des[idest[0]][idest[1]][idest[2]][ist] = factor*sou[ix][iy][iz][ist];
						 });
		
	} else {

		std::vector<long> point_list;

		for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < source.basis().sizes()[2]; iz++){
					
					auto ii = source.basis().to_symmetric_range(ix, iy, iz);
					auto idest = destination.basis().from_symmetric_range(ii);

					if(not destination.basis().local_contains(idest)) continue;

					point_list.push_back(source.basis().linear_index(ix, iy, iz));
					
				}
			}
		}

		gpu::array<long, 1> list(point_list.begin(), point_list.end());
		
		auto points = parallel::get_remote_points(source, list);
		
		long ip = 0;
		for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < source.basis().sizes()[2]; iz++){
	
					auto ii = source.basis().to_symmetric_range(ix, iy, iz);
					auto idest = destination.basis().from_symmetric_range(ii);

					if(not destination.basis().local_contains(idest)) continue;

					auto il0 = destination.basis().cubic_part(0).global_to_local(parallel::global_index(idest[0]));
					auto il1 = destination.basis().cubic_part(1).global_to_local(parallel::global_index(idest[1]));
					auto il2 = destination.basis().cubic_part(2).global_to_local(parallel::global_index(idest[2]));

					for(int ist = 0; ist < source.set_part().local_size(); ist++) destination.hypercubic()[il0][il1][il2][ist] = factor*points[ip][ist];
					ip++;
					
				}
			}
		}

	}
	
	return destination;			
}


//////////////////////////////////////////////////////////
		
template <class FieldType>
FieldType shrink(FieldType const & source, typename FieldType::basis_type const & new_basis, double const factor = 1.0) {

	CALI_CXX_MARK_FUNCTION;
	
	FieldType destination(new_basis);
			
	destination.fill(0.0);
		
	if(not new_basis.part().parallel()) {

		gpu::run(destination.basis().sizes()[2], destination.basis().sizes()[1], destination.basis().sizes()[0],
						 [sbas = source.basis().point_op(), dbas = destination.basis().point_op(), des = begin(destination.cubic()), sou = begin(source.cubic()), factor]
						 GPU_LAMBDA (auto iz, auto iy, auto ix){
							 
							 auto ii = dbas.to_symmetric_range(ix, iy, iz);
							 auto isource = sbas.from_symmetric_range(ii);
							 des[ix][iy][iz] = factor*sou[isource[0]][isource[1]][isource[2]];
						 });

	} else {

		gpu::array<long, 1> point_list(destination.basis().local_size());

		{
			long ip = 0;
			for(int ix = 0; ix < destination.basis().local_sizes()[0]; ix++){
				for(int iy = 0; iy < destination.basis().local_sizes()[1]; iy++){
					for(int iz = 0; iz < destination.basis().local_sizes()[2]; iz++){	

						auto ixg = destination.basis().cubic_part(0).local_to_global(ix);
						auto iyg = destination.basis().cubic_part(1).local_to_global(iy);
						auto izg = destination.basis().cubic_part(2).local_to_global(iz);						
						
						auto ii = destination.basis().to_symmetric_range(ixg, iyg, izg);
						auto isource = source.basis().from_symmetric_range(ii);
						
						point_list[ip] = source.basis().linear_index(isource[0], isource[1], isource[2]);
						ip++;
					}
				}
			}

			assert(ip == point_list.size());
		}

		auto points = parallel::get_remote_points(source, point_list);

		{
			long ip = 0;
			for(int ix = 0; ix < destination.basis().local_sizes()[0]; ix++){
				for(int iy = 0; iy < destination.basis().local_sizes()[1]; iy++){
					for(int iz = 0; iz < destination.basis().local_sizes()[2]; iz++){	
						
						destination.cubic()[ix][iy][iz] = factor*points[ip];
						ip++;
					}
				}
			}
		}
		
	}

	return destination;
}

//////////////////////////////////////////////////////////
		
template <class Type, class BasisType>
basis::field_set<BasisType, Type> shrink(basis::field_set<BasisType, Type> const & source, BasisType const & new_basis, double const factor = 1.0) {

	CALI_CXX_MARK_FUNCTION;
	
	basis::field_set<BasisType, Type> destination(new_basis, source.set_size(), source.full_comm());
			
	destination.fill(0.0);
	if(not new_basis.part().parallel()) {
			
		gpu::run(destination.basis().sizes()[2], destination.basis().sizes()[1], destination.basis().sizes()[0],
						 [nst = source.set_part().local_size(), sbas = source.basis().point_op(), dbas = destination.basis().point_op(), des = begin(destination.hypercubic()), sou = begin(source.hypercubic()), factor]
						 GPU_LAMBDA (auto iz, auto iy, auto ix){
							 
							 auto ii = dbas.to_symmetric_range(ix, iy, iz);
							 auto isource = sbas.from_symmetric_range(ii);
							 for(int ist = 0; ist < nst; ist++) des[ix][iy][iz][ist] = factor*sou[isource[0]][isource[1]][isource[2]][ist];
						 });
		
	} else {

		gpu::array<long, 1> point_list(destination.basis().local_size());
		
		{
			long ip = 0;
			for(int ix = 0; ix < destination.basis().local_sizes()[0]; ix++){
				for(int iy = 0; iy < destination.basis().local_sizes()[1]; iy++){
					for(int iz = 0; iz < destination.basis().local_sizes()[2]; iz++){	
						
						auto ixg = destination.basis().cubic_part(0).local_to_global(ix);
						auto iyg = destination.basis().cubic_part(1).local_to_global(iy);
						auto izg = destination.basis().cubic_part(2).local_to_global(iz);						
						
						auto ii = destination.basis().to_symmetric_range(ixg, iyg, izg);
						auto isource = source.basis().from_symmetric_range(ii);
						
						point_list[ip] = source.basis().linear_index(isource[0], isource[1], isource[2]);
						ip++;
					}
				}
			}

			assert(ip == point_list.size());
		}

		auto points = parallel::get_remote_points(source, point_list);

		{
			long ip = 0;
			for(int ix = 0; ix < destination.basis().local_sizes()[0]; ix++){
				for(int iy = 0; iy < destination.basis().local_sizes()[1]; iy++){
					for(int iz = 0; iz < destination.basis().local_sizes()[2]; iz++){	
						
						for(int ist = 0; ist < source.set_part().local_size(); ist++) destination.hypercubic()[ix][iy][iz][ist] = factor*points[ip][ist];
						ip++;
					}
				}
			}
		}

	}

	return destination;
}
		
//////////////////////////////////////////////////////////

template <class FieldType,
					typename std::enable_if<std::is_same<typename FieldType::element_type, complex>::value, int>::type = 0>
auto refine(FieldType const & source, typename basis::real_space const & new_basis){

	static_assert(std::is_same<typename FieldType::basis_type, basis::real_space>::value, "Only implemented for real space");
	static_assert(std::is_same<typename FieldType::element_type, complex>::value, "Only works for complex");
	
	assert(new_basis.size() == 8*source.basis().size()); //only a factor of 2 has been tested
	assert(not source.basis().part().parallel());
			
	basis::fourier_space new_fourier_basis(new_basis);
	auto destination_fourier = enlarge(operations::space::to_fourier(source), new_fourier_basis, 1.0/source.basis().size());
	return operations::space::to_real(destination_fourier,  /*normalize = */ false);
}

//////////////////////////////////////////////////////////

template <typename FieldType,
					typename std::enable_if<std::is_same<typename FieldType::element_type, double>::value, int>::type = 0>
auto refine(FieldType const & source, typename basis::real_space const & new_basis){

	auto complex_refine = refine(complex_field(source), new_basis);
	return real_field(complex_refine);
}

//////////////////////////////////////////////////////////
		
template <class FieldType,
					typename std::enable_if<std::is_same<typename FieldType::element_type, complex>::value, int>::type = 0>
auto coarsen(FieldType const & source, typename basis::real_space const & new_basis){

	assert(8*new_basis.size() == source.basis().size()); //only a factor of 2 has been tested		
	assert(not source.basis().part().parallel());
			
	basis::fourier_space new_fourier_basis(new_basis);
	auto destination_fourier = shrink(operations::space::to_fourier(source), new_fourier_basis, 1.0/source.basis().size());
	return operations::space::to_real(destination_fourier, /*normalize = */ false);
}

//////////////////////////////////////////////////////////

template <typename FieldType,
					typename std::enable_if<std::is_same<typename FieldType::element_type, double>::value, int>::type = 0>
auto coarsen(FieldType const & source, typename basis::real_space const & new_basis){

	auto complex_coarsen = coarsen(complex_field(source), new_basis);
	return real_field(complex_coarsen);
}

}

}
}
#endif

#ifdef INQ_OPERATIONS_TRANSFER_UNIT_TEST
#undef INQ_OPERATIONS_TRANSFER_UNIT_TEST

#include <operations/integral.hpp>

#include <catch2/catch_all.hpp>

namespace transfer_unit_test {
	
	auto outside(int ii, int ll){
		if(ii < 0.0){
			return ii < -(ll/2);
		} else {
			return ii >= (ll + 1)/2;
		}
	}
	
}

TEMPLATE_TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG, double, inq::complex) {

	using namespace transfer_unit_test;
	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

	
	parallel::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	auto set_comm = basis::set_subcomm(cart_comm);
	auto basis_comm = basis::basis_subcomm(cart_comm);

	parallel::communicator self_comm{boost::mpi3::environment::get_self_instance()};
	
	vector3<double> ll{6.66, 7.77, 9.99};
	
	systems::box box = systems::box::orthorhombic(ll[0]*1.0_b, ll[1]*1.0_b, ll[2]*1.0_b);
	auto spacing = 0.46320257;
	
	SECTION("Enlarge and shrink -- field"){
		
		basis::real_space grid(box, spacing, cart_comm);
		basis::field<basis::real_space, TestType> small(grid);
		
		CHECK(small.basis().rlength()[0] == Approx(ll[0]));
		CHECK(small.basis().rlength()[1] == Approx(ll[1]));
		CHECK(small.basis().rlength()[2] == Approx(ll[2]));
		
		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					auto rr = small.basis().point_op().rvector_cartesian(ix, iy, iz);
					small.cubic()[ix][iy][iz] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
				}
			}
		}

		auto large = operations::transfer::enlarge(small, grid.enlarge(2));

		CHECK(grid.part().parallel() == large.basis().part().parallel());
		CHECK(large.basis().rlength()[0] == Approx(2.0*ll[0]));
		CHECK(large.basis().rlength()[1] == Approx(2.0*ll[1]));
		CHECK(large.basis().rlength()[2] == Approx(2.0*ll[2]));

		long count_large = 0;
		long count_small = 0;
		for(int ix = 0; ix < large.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < large.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < large.basis().local_sizes()[2]; iz++){
				
					auto ixg = large.basis().cubic_part(0).local_to_global(ix);
					auto iyg = large.basis().cubic_part(1).local_to_global(iy);
					auto izg = large.basis().cubic_part(2).local_to_global(iz);						

					auto ii = large.basis().to_symmetric_range(ixg, iyg, izg);

					if(outside(ii[0], small.basis().sizes()[0]) or outside(ii[1], small.basis().sizes()[1]) or outside(ii[2], small.basis().sizes()[2])){
						CHECK(real(large.cubic()[ix][iy][iz]) == 0.0_a);
						CHECK(imag(large.cubic()[ix][iy][iz]) == 0.0_a);
						count_large++;
					} else {
						auto rr = large.basis().point_op().rvector_cartesian(ix, iy, iz);
						CHECK(real(large.cubic()[ix][iy][iz]) == Approx(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])));
						count_small++;
					}
				
				}
			}
		}

		cart_comm.all_reduce_n(&count_small, 1);
		cart_comm.all_reduce_n(&count_large, 1);

		CHECK(count_small == small.basis().size());
		CHECK(count_large > count_small);
		CHECK(count_large == large.basis().size() - count_small);

		auto small2 = operations::transfer::shrink(large, small.basis());

		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					CHECK(small2.cubic()[ix][iy][iz] == small.cubic()[ix][iy][iz]);
				}
			}
		}

	}

	SECTION("Enlarge and shrink -- field_set"){

		basis::real_space grid(box, spacing, basis_comm);
		basis::field_set<basis::real_space, TestType> small(grid, 5, cart_comm);
		
		CHECK(small.basis().rlength()[0] == Approx(ll[0]));
		CHECK(small.basis().rlength()[1] == Approx(ll[1]));
		CHECK(small.basis().rlength()[2] == Approx(ll[2]));
		
		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < small.set_part().local_size(); ist++){
							
						auto ixg = small.basis().cubic_part(0).local_to_global(ix);
						auto iyg = small.basis().cubic_part(1).local_to_global(iy);
						auto izg = small.basis().cubic_part(2).local_to_global(iz);						
						auto rr = small.basis().point_op().rvector_cartesian(ixg, iyg, izg);
						small.hypercubic()[ix][iy][iz][ist] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
					}
				}
			}
		}
		
		auto large = operations::transfer::enlarge(small, grid.enlarge(2));

		CHECK(large.basis().rlength()[0] == Approx(2.0*ll[0]));
		CHECK(large.basis().rlength()[1] == Approx(2.0*ll[1]));
		CHECK(large.basis().rlength()[2] == Approx(2.0*ll[2]));

		long count_large = 0;
		long count_small = 0;
		for(int ix = 0; ix < large.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < large.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < large.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < large.set_part().local_size(); ist++){
						
						auto ixg = large.basis().cubic_part(0).local_to_global(ix);
						auto iyg = large.basis().cubic_part(1).local_to_global(iy);
						auto izg = large.basis().cubic_part(2).local_to_global(iz);						
						
						auto ii = large.basis().to_symmetric_range(ixg, iyg, izg);
						
						if(outside(ii[0], small.basis().sizes()[0]) or outside(ii[1], small.basis().sizes()[1]) or outside(ii[2], small.basis().sizes()[2])){
							CHECK(real(large.hypercubic()[ix][iy][iz][ist]) == 0.0_a);
							CHECK(imag(large.hypercubic()[ix][iy][iz][ist]) == 0.0_a);							
							count_large++;
						} else {
							auto rr = large.basis().point_op().rvector_cartesian(ix, iy, iz);
							CHECK(real(large.hypercubic()[ix][iy][iz][ist]) == Approx(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])));
							count_small++;
						}
						
					}
				}
			}
		}
		
		cart_comm.all_reduce_n(&count_small, 1);
		cart_comm.all_reduce_n(&count_large, 1);

		CHECK(count_small == small.basis().size()*small.set_size());
		CHECK(count_large > count_small);
		CHECK(count_large == large.basis().size()*small.set_size() - count_small);
	
		auto small2 = operations::transfer::shrink(large, small.basis());
	
		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < small.set_part().local_size(); ist++){
						CHECK(small2.hypercubic()[ix][iy][iz][ist] == small.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
	}

	SECTION("Mesh refinement -- field"){

		basis::real_space grid(box, spacing, self_comm);
		basis::field<basis::real_space, TestType> coarse(grid); 

		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					
					auto ixg = coarse.basis().cubic_part(0).local_to_global(ix);
					auto iyg = coarse.basis().cubic_part(1).local_to_global(iy);
					auto izg = coarse.basis().cubic_part(2).local_to_global(iz);						
					auto rr = coarse.basis().point_op().rvector_cartesian(ixg, iyg, izg);
					coarse.cubic()[ix][iy][iz] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
				}
			}
		}

		auto fine = operations::transfer::refine(coarse, grid.refine(2));

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		CHECK(real(operations::integral(coarse)) == 109.2100704359_a);
		CHECK(fabs(imag(operations::integral(coarse))) < 1e-14);
		CHECK(real(operations::integral(fine)) == 109.2100704359_a);
		CHECK(fabs(imag(operations::integral(fine))) < 1e-14);

			
		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					
					auto ixg = fine.basis().cubic_part(0).local_to_global(ix);
					auto iyg = fine.basis().cubic_part(1).local_to_global(iy);
					auto izg = fine.basis().cubic_part(2).local_to_global(iz);						
					auto rr = fine.basis().point_op().rvector_cartesian(ixg, iyg, izg);
					CHECK(fabs(real(fine.cubic()[ix][iy][iz])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
					CHECK(fabs(imag(fine.cubic()[ix][iy][iz])) < 5e-3);
				}
			}
		}
		
	}
		
	SECTION("Mesh refinement -- field_set"){

		basis::real_space grid(box, spacing, self_comm);
		basis::field_set<basis::real_space, TestType> coarse(grid, 5); 

		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < coarse.set_part().local_size(); ist++){
						
						auto ixg = coarse.basis().cubic_part(0).local_to_global(ix);
						auto iyg = coarse.basis().cubic_part(1).local_to_global(iy);
						auto izg = coarse.basis().cubic_part(2).local_to_global(iz);						
						auto rr = coarse.basis().point_op().rvector_cartesian(ixg, iyg, izg);
						coarse.hypercubic()[ix][iy][iz][ist] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);

					}
				}
			}
		}

		auto fine = operations::transfer::refine(coarse, grid.refine(2));

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < fine.set_part().local_size(); ist++){
							
						auto ixg = fine.basis().cubic_part(0).local_to_global(ix);
						auto iyg = fine.basis().cubic_part(1).local_to_global(iy);
						auto izg = fine.basis().cubic_part(2).local_to_global(iz);						
						auto rr = fine.basis().point_op().rvector_cartesian(ixg, iyg, izg);
						CHECK(fabs(real(fine.hypercubic()[ix][iy][iz][ist])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
						CHECK(fabs(imag(fine.hypercubic()[ix][iy][iz][ist])) < 5e-3);
						
					}
				}
			}
		}
		
	}
	
	SECTION("Mesh coarsening -- field"){

		basis::real_space grid(box, spacing, self_comm);
		
		auto fine_grid = grid.refine(2);
		
		basis::field<basis::real_space, TestType> fine(fine_grid); 

		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					
					auto ixg = fine.basis().cubic_part(0).local_to_global(ix);
					auto iyg = fine.basis().cubic_part(1).local_to_global(iy);
					auto izg = fine.basis().cubic_part(2).local_to_global(iz);						
					auto rr = fine.basis().point_op().rvector_cartesian(ixg, iyg, izg);
					fine.cubic()[ix][iy][iz] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
				}
			}
		}

		auto coarse = operations::transfer::coarsen(fine, grid);

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		CHECK(real(operations::integral(coarse)) == 109.3027787767_a);
		CHECK(fabs(imag(operations::integral(coarse))) < 1e-14);
		CHECK(real(operations::integral(fine)) == 109.3027787767_a);
		CHECK(fabs(imag(operations::integral(fine))) < 1e-14);

			
		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					
					auto ixg = coarse.basis().cubic_part(0).local_to_global(ix);
					auto iyg = coarse.basis().cubic_part(1).local_to_global(iy);
					auto izg = coarse.basis().cubic_part(2).local_to_global(iz);						
					auto rr = coarse.basis().point_op().rvector_cartesian(ixg, iyg, izg);
					CHECK(fabs(real(coarse.cubic()[ix][iy][iz])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
					CHECK(fabs(imag(coarse.cubic()[ix][iy][iz])) < 5e-3);
				}
			}
		}
		
	}

	SECTION("Mesh coarsening -- field_set"){
			
		basis::real_space grid(box, spacing, self_comm);
		auto fine_grid = grid.refine(2);
		
		basis::field_set<basis::real_space, TestType> fine(fine_grid, 5);

		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < fine.set_part().local_size(); ist++){
						
						auto ixg = fine.basis().cubic_part(0).local_to_global(ix);
						auto iyg = fine.basis().cubic_part(1).local_to_global(iy);
						auto izg = fine.basis().cubic_part(2).local_to_global(iz);						
						auto rr = fine.basis().point_op().rvector_cartesian(ixg, iyg, izg);
						fine.hypercubic()[ix][iy][iz][ist] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
					}
				}
			}
		}

		auto coarse = operations::transfer::coarsen(fine, grid);

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < fine.set_part().local_size(); ist++){
						
						auto ixg = coarse.basis().cubic_part(0).local_to_global(ix);
						auto iyg = coarse.basis().cubic_part(1).local_to_global(iy);
						auto izg = coarse.basis().cubic_part(2).local_to_global(iz);						
						auto rr = coarse.basis().point_op().rvector_cartesian(ixg, iyg, izg);
						CHECK(fabs(real(coarse.hypercubic()[ix][iy][iz][ist])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
						CHECK(fabs(imag(coarse.hypercubic()[ix][iy][iz][ist])) < 5e-3);
					}
				}
			}
		}
		
	}

}
#endif


