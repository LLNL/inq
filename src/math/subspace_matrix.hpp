/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__SUBSPACE_MATRIX
#define INQ__MATH__SUBSPACE_MATRIX

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa.

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

#include <math/array.hpp>
#include <gpu/run.hpp>
#include <utils/partition.hpp>

#include <mpi3/cartesian_communicator.hpp>

namespace inq {
namespace math {

template <class Type>
class subspace_matrix {

public:
  
  using array_type = math::array<Type, 2>;
  
  subspace_matrix(boost::mpi3::cartesian_communicator<2> & comm, long size):
    comm_(comm),
    array_({size, size}),
		part_(size, comm.axis(0)){
  }

  subspace_matrix(boost::mpi3::cartesian_communicator<2> & comm, long size, Type const & ival):
    comm_(comm),
    array_({size, size}, ival),
		part_(size, comm.axis(0)){
  }
  
  subspace_matrix(boost::mpi3::cartesian_communicator<2> & comm, array_type && mat):
    comm_(comm),
    array_(std::move(mat)),
		part_(array_.size(), comm.axis(0)){
  }

  auto size() const {
    return array_.size();
  }
  
  auto & array() const {
    return array_;
  }

  auto & array() {
    return array_;
  }

  math::array<Type, 1> diagonal() const {
    math::array<Type, 1> diag(part_.local_size());
    gpu::run(part_.local_size(), [dia = begin(diag), arr = begin(array_), pa = part_] GPU_LAMBDA (auto ii){
			auto iig = pa.start() + ii;
			dia[ii] = arr[iig][iig];
		});
    return diag;
  }

  auto comm() const {
    return comm_;
  }
  
private:

  mutable boost::mpi3::cartesian_communicator<2> comm_;
  array_type array_;
	utils::partition part_;
  
};

}
}

#ifdef INQ_MATH_SUBSPACE_MATRIX_UNIT_TEST
#undef INQ_MATH_SUBSPACE_MATRIX_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("math::subspace_matrix", "[math::subspace_matrix]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  auto comm = boost::mpi3::environment::get_world_instance();
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {});
  
  math::array<double, 2> matrix({2, 2});
  
  matrix[0][0] = 4.0;
  matrix[0][1] = 0.0;
  matrix[1][0] = 0.0;
  matrix[1][1] = 2.0;

  math::subspace_matrix mm(cart_comm, std::move(matrix));

  CHECK(mm.array()[0][0] == 4.0);
  CHECK(mm.array()[0][1] == 0.0);
  CHECK(mm.array()[1][0] == 0.0);
  CHECK(mm.array()[1][1] == 2.0);

  auto dd = mm.diagonal();

  CHECK(dd[0] == 4.0);
  CHECK(dd[1] == 2.0);
  
}

#endif
#endif
