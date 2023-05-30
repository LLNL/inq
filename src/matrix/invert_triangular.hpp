/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATRIX__INVERT_TRIANGULAR
#define INQ__MATRIX__INVERT_TRIANGULAR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/array.hpp>
#include <matrix/gather_scatter.hpp>
#include <solvers/invert_triangular.hpp>

namespace inq {
namespace matrix {

template <typename DistributedMatrix>
void invert_triangular(DistributedMatrix & matrix) {

  assert(matrix.sizex() == matrix.sizey());

  using type = typename DistributedMatrix::element_type;
  
  auto full_matrix = matrix::gather(matrix, /* root = */ 0);

	static_assert(std::is_same_v<type, double> or std::is_same_v<type, complex>, "invert_triangular is only implemented for double and complex");

  if(matrix.comm().root()){
    int nn = matrix.sizex();
    
    gpu::array<type, 2> inverse({nn, nn});
    
    gpu::run(nn, nn, [inv = begin(inverse)] GPU_LAMBDA (auto jj, auto ii){
      inv[ii][jj] = (ii == jj) ? 1.0 : 0.0;
    });
    
    namespace blas = boost::multi::blas;
    blas::trsm(blas::side::right, blas::filling::upper, 1.0, blas::H(full_matrix), inverse);
    
    gpu::run(nn, nn, [mat = begin(full_matrix), inv = begin(inverse)] GPU_LAMBDA (auto ii, auto jj){
      mat[ii][jj] = conj(inv[jj][ii]);
    });
  }
  
  matrix::scatter(full_matrix, matrix, /* root = */ 0);

}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_MATRIX_INVERT_TRIANGULAR_UNIT_TEST
#undef INQ_MATRIX_INVERT_TRIANGULAR_UNIT_TEST

#include <math/complex.hpp>
#include <catch2/catch_all.hpp>

using namespace inq;
using namespace Catch::literals;

TEMPLATE_TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG, double, complex) {

  parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	inq::parallel::cartesian_communicator<2> cart_comm(comm, {});

	SECTION("2x2"){
	
    gpu::array<complex, 2> array({2, 2});
  		
		array[0][0] = 4.0;
		array[1][0] = -1.0;
		array[1][1] = 2.0;

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		matrix::invert_triangular(matrix);

    array = matrix::all_gather(matrix);

		CHECK(array[0][0] == 0.25);
		CHECK(array[0][1] == 0.0);
		CHECK(array[1][0] == 0.125);
		CHECK(array[1][1] == 0.5);
  }

	SECTION("NxN"){
	
		auto nn = 15;

    gpu::array<complex, 2> array({nn, nn});
		
		for(int ii = 0; ii < nn; ii++){
			for(int jj = 0; jj < nn; jj++){
				if(ii < jj) {
					array[ii][jj] = 0.0;
				} else {
					array[ii][jj] = cos(ii)*(jj + 0.1) + sin(jj - 0.25)*(ii + 1.0);
				}
				if constexpr (std::is_same_v<TestType, complex>) array[ii][jj]*= exp(complex(0.0, (ii + 1.0)*(cos(jj) - 0.25)));
			}
		}

    matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
      
		matrix::invert_triangular(matrix);

    auto inverse = matrix::all_gather(matrix);

		auto mul = +boost::multi::blas::gemm(1.0, inverse, array);
    
		for(int ii = 0; ii < nn; ii++){
			for(int jj = 0; jj < nn; jj++){
				if(ii == jj) {
					CHECK(real(mul[ii][jj]) == 1.0_a);
					CHECK(imag(mul[ii][jj]) < 1e-12);					
				} else {
					CHECK(fabs(mul[ii][jj]) < 1e-12);
				}
			}
		}
			
  }
}

#endif
