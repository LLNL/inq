/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATRIX__CHOLESKY
#define INQ__MATRIX__CHOLESKY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <matrix/gather_scatter.hpp>
#include <solvers/cholesky.hpp>

namespace inq {
namespace matrix {

template <typename DistributedMatrix>
void cholesky(DistributedMatrix & matrix) {

  assert(matrix.sizex() == matrix.sizey());
  
  auto full_matrix = matrix::gather(matrix, /* root = */ 0);
  if(matrix.comm().root()) solvers::cholesky(full_matrix);
  matrix::scatter(full_matrix, matrix, /* root = */ 0);

}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_MATRIX_CHOLESKY_UNIT_TEST
#undef INQ_MATRIX_CHOLESKY_UNIT_TEST

#include <gpu/array.hpp>
#include <matrix/distributed.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
  using namespace Catch::literals;
  
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	parallel::cartesian_communicator<2> cart_comm(comm, {});


  SECTION("Complex 2x2"){
    
		using namespace inq;
    using namespace Catch::literals;
		
		gpu::array<complex, 2> array({2, 2});
		
		array[0][0] = 6432.12;
		array[0][1] = 4502.48;
		array[1][0] = 4502.48;
		array[1][1] = 3151.74;

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		matrix::cholesky(matrix);

    array = matrix::all_gather(matrix);
	 
		CHECK(real(array[0][0]) == 80.2005_a);
		CHECK(real(array[0][1]) == 0.0_a);
		CHECK(real(array[1][0]) == 56.1402992511_a);
		CHECK(real(array[1][1]) == 0.0824620974_a);    
  }

}
#endif
