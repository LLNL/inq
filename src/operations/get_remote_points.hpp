/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__GET_REMOTE_POINTS
#define INQ__UTILS__GET_REMOTE_POINTS

/*
 Copyright (C) 2019 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <basis/field.hpp>

#include <mpi3/communicator.hpp>

namespace inq {
namespace operations {

template <class BasisType, class ElementType, class ArrayType>
auto get_remote_points(basis::field<BasisType, ElementType> const & source, ArrayType const & point_list){

	CALI_CXX_MARK_FUNCTION;

	struct point_position {
		long point;
		long position;
	};

	auto const num_proc = source.basis().comm().size();

	// create a list per processor of the points we need
	std::vector<std::vector<point_position>> points_needed(num_proc);
		
	for(long ilist = 0; ilist < point_list.size(); ilist++){
		auto src_proc = source.basis().part().location(point_list[ilist]);
		assert(src_proc >= 0 and src_proc < num_proc); 
		points_needed[src_proc].push_back(point_position{point_list[ilist], ilist});
	}

	// find out how many points each processor requests from us
	math::array<int, 1> list_sizes_needed(num_proc);
	math::array<int, 1> list_sizes_requested(num_proc);	
	
	int total_needed = 0;
	for(int iproc = 0; iproc < num_proc; iproc++){
		list_sizes_needed[iproc] = points_needed[iproc].size();
		total_needed += list_sizes_needed[iproc];
	}

	assert(total_needed == point_list.size());

	MPI_Alltoall(list_sizes_needed.data(), 1, MPI_INT, list_sizes_requested.data(), 1, MPI_INT, source.basis().comm().get());

	// get the list of points each processor requests from us and send a list of the points we need
	math::array<int, 1> list_needed(total_needed);
	math::array<int, 1> list_displs_needed(num_proc);
	math::array<int, 1> list_displs_requested(num_proc);
	
	int total_requested = 0;

	source.basis().comm().barrier();
	
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

	source.basis().comm().barrier();
	
	for(int iproc = 0; iproc < num_proc; iproc++){

		source.basis().comm().barrier();
		
		for(long ip = 0; ip < list_sizes_needed[iproc]; ip++) {
			list_needed[list_displs_needed[iproc] + ip] = points_needed[iproc][ip].point;
		}
	}

	math::array<int, 1> list_points_requested(total_requested);
	
	MPI_Alltoallv(list_needed.data(), list_sizes_needed.data(), list_displs_needed.data(), MPI_INT, list_points_requested.data(), list_sizes_requested.data(), list_displs_requested.data(), MPI_INT, source.basis().comm().get());

	// send the value of the request points
	math::array<ElementType, 1> value_points_requested(total_requested);
	math::array<ElementType, 1> value_points_needed(total_needed);
	
	for(long ip = 0; ip < total_requested; ip++){
		auto iplocal = source.basis().part().global_to_local(utils::global_index(list_points_requested[ip]));
		assert(iplocal < source.basis().size());
		value_points_requested[ip] = source.linear()[iplocal];
	}

	auto mpi_type = boost::mpi3::detail::basic_datatype<ElementType>();
	
	MPI_Alltoallv(value_points_requested.data(), list_sizes_requested.data(), list_displs_requested.data(), mpi_type,
								value_points_needed.data(), list_sizes_needed.data(), list_displs_needed.data(), mpi_type, source.basis().comm().get());

	// Finally copy the values to the return array in the proper order
	math::array<ElementType, 1> remote_points(point_list.size());

	long ip = 0;
	for(int iproc = 0; iproc < num_proc; iproc++){
		for(long jp = 0; jp < long(points_needed[iproc].size()); jp++) {
			remote_points[points_needed[iproc][jp].position] = value_points_needed[ip];
			ip++;
		}
	}

	assert(ip == point_list.size());
	
  return remote_points;
}

}
}

#ifdef INQ_OPERATIONS_GET_REMOTE_POINTS_UNIT_TEST
#undef INQ_OPERATIONS_GET_REMOTE_POINTS_UNIT_TEST

#include <math/vector3.hpp>

#include <catch2/catch.hpp>
#include <mpi3/cartesian_communicator.hpp>

TEST_CASE("Class operations::get_remote_points", "[operations::get_remote_points]"){

	using namespace inq;
	using namespace Catch::literals;
	using math::vector3;

  boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

	double lx = 13.3;
	double ly = 6.55;
	double lz = 8.02;
 	ions::UnitCell cell(vector3<double>(lx, 0.0, 0.0), vector3<double>(0.0, ly, 0.0), vector3<double>(0.0, 0.0, lz));
  basis::real_space rs(cell, input::basis::cutoff_energy(20.0), cart_comm);

  basis::field<basis::real_space, complex> test_field(rs);

  for(long ip = 0; ip < rs.local_size(); ip++){
    auto ipg = rs.part().local_to_global(ip);
    test_field.linear()[ip] = complex(ipg.value(), 0.1*ipg.value());
  }

	srand48(500 + 34895783*cart_comm.rank());
	
	//	long const npoints = 1;//drand48()*rs.size();
	long const npoints = drand48()*rs.size();

	math::array<long, 1> list(npoints);
	
	for(long ip = 0; ip < npoints; ip++){
		list[ip] = drand48()*(rs.size() - 1);
		assert(list[ip] < rs.size());
	}

	auto remote_points = operations::get_remote_points(test_field, list);

	for(long ip = 0; ip < npoints; ip++){
		CHECK(double(list[ip]) == real(remote_points[ip]));
		CHECK(0.1*double(list[ip]) == imag(remote_points[ip]));
	}

	
}


#endif
#endif


