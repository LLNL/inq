/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__INVERT_TRIANGULAR
#define INQ__SOLVERS__INVERT_TRIANGULAR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <FC.h>

#include <tuple> //std::get
#include <cassert>

#include <gpu/array.hpp>
#include <math/complex.hpp>
#include <math/subspace_matrix.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace solvers {

template <typename Type>
void invert_triangular(math::subspace_matrix<Type> & matrix){
  CALI_CXX_MARK_SCOPE("invert_triangular_double");

	static_assert(std::is_same_v<Type, double> or std::is_same_v<Type, complex>, "invert_triangular is only implemented for double and complex");
	
	int nn = std::get<0>(sizes(matrix.array()));

	gpu::array<Type, 2> inverse({nn, nn});

	gpu::run(nn, nn, [inv = begin(inverse)] GPU_LAMBDA (auto jj, auto ii){
		inv[ii][jj] = (ii == jj) ? 1.0 : 0.0;
	});

	namespace blas = boost::multi::blas;
	blas::trsm(blas::side::right, blas::filling::upper, 1.0, blas::H(matrix.array()), inverse);

	gpu::run(nn, nn, [mat = begin(matrix.array()), inv = begin(inverse)] GPU_LAMBDA (auto ii, auto jj){
		mat[ii][jj] = conj(inv[jj][ii]);
	});

}

}
}
#endif

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_INVERT_TRIANGULAR_UNIT_TEST
#undef INQ_SOLVERS_INVERT_TRIANGULAR_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <mpi3/environment.hpp>

#include <gpu/array.hpp>

TEMPLATE_TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG, double, complex) {

  parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	inq::parallel::cartesian_communicator<2> cart_comm(comm, {});

	SECTION("2x2"){
	
		using namespace inq;
		using namespace Catch::literals;

    math::subspace_matrix<TestType> matrix(cart_comm, 2);
		
		matrix.array()[0][0] = 4.0;
		matrix.array()[1][0] = -1.0;
		matrix.array()[1][1] = 2.0;

		solvers::invert_triangular(matrix);
    
    CHECK(matrix.array()[0][0] == 0.25);
		CHECK(matrix.array()[0][1] == 0.0);
    CHECK(matrix.array()[1][0] == 0.125);
    CHECK(matrix.array()[1][1] == 0.5);
  }

	SECTION("NxN"){
	
		using namespace inq;
		using namespace Catch::literals;

		auto nn = 15;

    math::subspace_matrix<TestType> matrix(cart_comm, nn);
		
		for(int ii = 0; ii < nn; ii++){
			for(int jj = 0; jj < nn; jj++){
				if(ii < jj) {
					matrix.array()[ii][jj] = 0.0;
				} else {
					matrix.array()[ii][jj] = cos(ii)*(jj + 0.1) + sin(jj - 0.25)*(ii + 1.0);
				}
				if constexpr (std::is_same_v<TestType, complex>) matrix.array()[ii][jj]*= exp(complex(0.0, (ii + 1.0)*(cos(jj) - 0.25)));
			}
		}

		auto orig = matrix;
		
		solvers::invert_triangular(matrix);
		
		auto mul = +boost::multi::blas::gemm(1.0, matrix.array(), orig.array());

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
