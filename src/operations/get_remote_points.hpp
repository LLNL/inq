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

	math::array<ElementType, 1> remote_points(point_list.size());
	
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

  for(long ip = 0; ip < rs.size(); ip++){
    auto ipg = rs.part().local_to_global(ip);
    test_field.linear()[ip] = complex(ipg.value(), 0.1*ipg.value());
  }

	srand48(500 + 34895783*cart_comm.rank());
	
	long const npoints = drand48()*rs.size();

	std::cout << "NPOINTS " << npoints << std::endl;
	
	std::vector<long> list(npoints);
	
	for(long ip = 0; ip < npoints; ip++){
		list[ip] = drand48()*(rs.size() - 1);
		assert(list[ip] < rs.size());
	}

	auto remote_points = operations::get_remote_points(test_field, list);

	for(long ip = 0; ip < npoints; ip++){
		CHECK(double(list[ip]) == real(remote_points[ip]));
		CHECK(double(list[ip]) == imag(0.1*remote_points[ip]));
	}

	
}


#endif
#endif


