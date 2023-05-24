/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__LEAST_SQUARES
#define INQ__SOLVERS__LEAST_SQUARES

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/array.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <FC.h>

#include <stdexcept>
#include <tuple> //for get

#define dgelsd FC_GLOBAL(dgelsd, DGELSD) 
extern "C" void dgelsd(const int & m, const int & n, const int & nrhs, double * a, const int & lda, double * b, const int & ldb, double * s, const double & rcond, int & rank, double * work, const int & lwork, int * iwork, int & info);

namespace inq {
namespace solvers {

template <class matrix_type, class vector_type>
void least_squares(matrix_type && matrix, vector_type & rhs){

	CALI_CXX_MARK_FUNCTION;
	
	int mm = std::get<0>(sizes(matrix));
	int nn = std::get<1>(sizes(matrix));

	gpu::array<double, 1> ss(mm);

	int rank, info;
	double dwork;
	int iwork_size;
	
	dgelsd(mm, nn, 1, raw_pointer_cast(matrix.data_elements()), mm, raw_pointer_cast(rhs.data_elements()), mm, raw_pointer_cast(ss.base()), -1.0, rank, &dwork, -1, &iwork_size, info);

	if(info != 0) std::runtime_error("inq error: dgelss in least_squares failed with info = " + std::to_string(info));
	
	auto work = (double *) malloc(int(dwork)*sizeof(double));
	auto iwork = (int *) malloc(iwork_size*sizeof(int));

	dgelsd(mm, nn, 1, raw_pointer_cast(matrix.data_elements()), mm, raw_pointer_cast(rhs.data_elements()), mm, raw_pointer_cast(ss.data_elements()), -1.0, rank, work, int(dwork), iwork, info);

	if(info != 0) std::runtime_error("inq error: dgelss in least_squares failed with info = " + std::to_string(info));
		
	free(work);
	free(iwork);
	
}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_LEAST_SQUARES_UNIT_TEST
#undef INQ_SOLVERS_LEAST_SQUARES_UNIT_TEST

#include <catch2/catch_all.hpp>

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
		
		solvers::least_squares(matrix, vector);

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
		
		solvers::least_squares(matrix, vector);

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
		
		solvers::least_squares(matrix, vector);

		CHECK(vector[0] == -30.882245351_a);
		CHECK(vector[1] ==  44.1177546523_a);
		
  }

}
#endif
