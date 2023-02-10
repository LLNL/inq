/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__LINEAR
#define INQ__SOLVERS__LINEAR

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
#include <utils/raw_pointer_cast.hpp>

#define dpotrf FC_GLOBAL(dpotrf, DPOTRF) 
extern "C" void dpotrf(const char * uplo, const int * n, double * a, const int * lda, int * info);

#define dtrsv FC_GLOBAL(dtrsv, DTRSV) 

namespace inq {
namespace solvers {

  /* This function calculates the inverse of a matrix by
     diagonalization, it is slow but accurate. */
  
template <class matrix_type, class vector_type>
void linear_symmetric(matrix_type && matrix, vector_type & vector){

	// the matrix must be square
	assert(std::get<0>(sizes(matrix)) == std::get<1>(sizes(matrix)));

	int nn = std::get<0>(sizes(matrix));
    
	int info;
	dpotrf("U", &nn, raw_pointer_cast(matrix.data_elements()), &nn, &info);
		
	const int one = 1;
	dtrsv('U', 'T', 'N', nn, raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(vector.data_elements()), one);
	dtrsv('U', 'N', 'N', nn, raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(vector.data_elements()), one);
	
}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_LINEAR_UNIT_TEST
#undef INQ_SOLVERS_LINEAR_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <math/array.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	SECTION("Diagonal real 2x2"){
	
		using namespace inq;
		using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 4.0;
		matrix[0][1] = 0.0;
		matrix[1][0] = 0.0;
		matrix[1][1] = 2.0;

		math::array<double, 1> vector = {0.0, 1.0, };
		
		solvers::linear_symmetric(matrix, vector);

		CHECK(vector[0] == 0.0_a);
		CHECK(vector[1] == 0.5_a);

  }
	
	SECTION("Symmetric real 2x2 -- 1"){
	
		using namespace inq;
	using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 0.89653;
		matrix[0][1] = 0.41072;
		matrix[1][0] = 0.41072;
		matrix[1][1] = 0.69479;

		math::array<double, 1> vector = {0.21563, 0.40103, };
		
		solvers::linear_symmetric(matrix, vector);

		CHECK(vector[0] == -0.0327899021_a);
		CHECK(vector[1] == 0.596579_a);
		
  }
	
	SECTION("Symmetric real 2x2 -- 2"){
	
		using namespace inq;
	using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 6432.12;
		matrix[0][1] = 4502.48;
		matrix[1][0] = 4502.48;
		matrix[1][1] = 3151.74;

		math::array<double, 1> vector({1.0, 1.0});
		
		solvers::linear_symmetric(matrix, vector);

		CHECK(vector[0] == -30.882245351_a);
		CHECK(vector[1] ==  44.1177546523_a);
		
  }
}
#endif
