/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__SUBSPACE_MATRIX
#define INQ__MATH__SUBSPACE_MATRIX

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/array.hpp>
#include <gpu/run.hpp>
#include <parallel/partition.hpp>

#include <parallel/communicator.hpp>

namespace inq {
namespace math {

template <class Type>
class subspace_matrix {

public:
  
  using array_type = gpu::array<Type, 2>;
  
  subspace_matrix(parallel::cartesian_communicator<2> & comm, long size):
    comm_(comm),
    array_({size, size}),
		part_(size, comm.axis(1)){
  }

  subspace_matrix(parallel::cartesian_communicator<2> & comm, long size, Type const & ival):
    comm_(comm),
    array_({size, size}, ival),
		part_(size, comm.axis(1)){
  }
  
  subspace_matrix(parallel::cartesian_communicator<2> & comm, array_type && mat):
    comm_(comm),
    array_(std::move(mat)),
		part_(array_.size(), comm.axis(1)){
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

  gpu::array<Type, 1> diagonal() const {
    gpu::array<Type, 1> diag(part_.local_size());
    gpu::run(part_.local_size(), [dia = begin(diag), arr = begin(array_), pa = part_] GPU_LAMBDA (auto ii){
			auto iig = pa.start() + ii;
			dia[ii] = arr[iig][iig];
		});
    return diag;
  }

  auto comm() const {
    return comm_;
  }

	auto & part() const {
		return part_;
	}
	
private:

  mutable parallel::cartesian_communicator<2> comm_;
  array_type array_;
	parallel::partition part_;
  
};

}
}
#endif

#ifdef INQ_MATH_SUBSPACE_MATRIX_UNIT_TEST
#undef INQ_MATH_SUBSPACE_MATRIX_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	parallel::cartesian_communicator<2> cart_comm(comm, {});
  
  gpu::array<double, 2> matrix({2, 2});
  
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

  if(mm.part().contains(0)) CHECK(dd[mm.part().global_to_local(parallel::global_index(0))] == 4.0);
  if(mm.part().contains(1)) CHECK(dd[mm.part().global_to_local(parallel::global_index(1))] == 2.0);	
  
}
#endif
