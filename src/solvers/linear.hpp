/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__LINEAR
#define INQ__SOLVERS__LINEAR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <FC.h>

#include <tuple> //std::get
#include <cassert>
#include <gpu/array.hpp>
#include <utils/raw_pointer_cast.hpp>

#define dpotrf FC_GLOBAL(dpotrf, DPOTRF) 
extern "C" void dpotrf(const char * uplo, const int * n, double * a, const int * lda, int * info);

#define dtrsv FC_GLOBAL(dtrsv, DTRSV) 

namespace inq {
namespace solvers {

template <class matrix_type, class vector_type>
void linear_symmetric(matrix_type && matrix, vector_type & vector){

	using std::get;

	assert(get<0>(sizes(matrix)) == get<1>(sizes(matrix)));

	int nn = get<0>(sizes(matrix));
    
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

#include <gpu/array.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	SECTION("Diagonal real 2x2"){
	
		using namespace inq;
		using namespace Catch::literals;
		
		gpu::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 4.0;
		matrix[0][1] = 0.0;
		matrix[1][0] = 0.0;
		matrix[1][1] = 2.0;

		gpu::array<double, 1> vector = {0.0, 1.0, };
		
		solvers::linear_symmetric(matrix, vector);

		CHECK(vector[0] == 0.0_a);
		CHECK(vector[1] == 0.5_a);

  }
	
	SECTION("Symmetric real 2x2 -- 1"){
	
		using namespace inq;
	using namespace Catch::literals;
		
		gpu::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 0.89653;
		matrix[0][1] = 0.41072;
		matrix[1][0] = 0.41072;
		matrix[1][1] = 0.69479;

		gpu::array<double, 1> vector = {0.21563, 0.40103, };
		
		solvers::linear_symmetric(matrix, vector);

		CHECK(vector[0] == -0.0327899021_a);
		CHECK(vector[1] == 0.596579_a);
		
  }
	
	SECTION("Symmetric real 2x2 -- 2"){
	
		using namespace inq;
	using namespace Catch::literals;
		
		gpu::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 6432.12;
		matrix[0][1] = 4502.48;
		matrix[1][0] = 4502.48;
		matrix[1][1] = 3151.74;

		gpu::array<double, 1> vector({1.0, 1.0});
		
		solvers::linear_symmetric(matrix, vector);

		CHECK(vector[0] == -30.882245351_a);
		CHECK(vector[1] ==  44.1177546523_a);
		
  }
}
#endif
