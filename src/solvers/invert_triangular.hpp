/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__INVERT_TRIANGULAR
#define INQ__SOLVERS__INVERT_TRIANGULAR

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

#include <FC.h>

#include <tuple> //std::get
#include <cassert>

#include <math/array.hpp>
#include <math/complex.hpp>
#include <math/subspace_matrix.hpp>
#include <utils/raw_pointer_cast.hpp>

#define dtrtri FC_GLOBAL(dtrtri, DTRTRI) 
extern "C" void  dtrtri(const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);

#define ztrtri FC_GLOBAL(ztrtri, ZTRTRI) 
extern "C" void  ztrtri(const char * uplo, const char * diag, const int * n, inq::complex * a, const int * lda, int * info);

namespace inq {
namespace solvers {

template <typename Type>
void invert_triangular(math::subspace_matrix<Type> & matrix){
  CALI_CXX_MARK_SCOPE("invert_triangular_double");

	static_assert(std::is_same_v<Type, double> or std::is_same_v<Type, complex>, "invert_triangular is only implemented for double and complex");
	
	int nn = std::get<0>(sizes(matrix.array()));

	auto matrix_data = raw_pointer_cast(matrix.array().data_elements());
	int info;
	
	if constexpr (std::is_same_v<Type, double>){
		dtrtri("U", "N", &nn, matrix_data, &nn, &info);
	} else {
		ztrtri("U", "N", &nn, matrix_data, &nn, &info);
	}

	gpu::run(nn, nn, [mat = begin(matrix.array())] GPU_LAMBDA (auto ii, auto jj){
		if(ii < jj) mat[ii][jj] = 0.0;
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

#include <math/array.hpp>

TEMPLATE_TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG, double, complex) {

  auto comm = boost::mpi3::environment::get_world_instance();
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
