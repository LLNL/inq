/* -*- indent-tabs-mode: t -*- */

#ifndef SOLVERS__LEAST_SQUARES
#define SOLVERS__LEAST_SQUARES

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <config.h>
#include <cstdlib>

#define dgelss FC_FUNC(dgelss, DGELSS) 
extern "C" void dgelss(const int & m, const int & n, const int & nrhs, double * a, const int & lda, double * b, const int & ldb, double * s, const double & rcond, int & rank, double * work, const int & lwork, int & info);

namespace solvers {

  template <class matrix_type, class vector_type>
  void least_squares(matrix_type && matrix, vector_type & rhs){

    int mm = std::get<0>(sizes(matrix));
		int nn = std::get<1>(sizes(matrix));

		math::array<double, 1> ss(mm);

		int rank, info;
		double dwork;
		
		//DATAOPERATIONS RAWLAPACK dgelss
		dgelss(mm, nn, 1, matrix.data(), mm, rhs.data(), mm, ss.data(), -1.0, rank, &dwork, -1, info);

		assert(info == 0);
		
		auto work = (double *) malloc(int(dwork)*sizeof(double));

		dgelss(mm, nn, 1, matrix.data(), mm, rhs.data(), mm, ss.data(), -1.0, rank, work, int(dwork), info);

		assert(info == 0);
		
		free(work);
		
  }

}

///////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

#include <math/array.hpp>

TEST_CASE("function solvers::least_squares", "[solvers::least_squares]") {

	SECTION("Diagonal real 2x2"){
	
		using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 4.0;
		matrix[0][1] = 0.0;
		matrix[1][0] = 0.0;
		matrix[1][1] = 2.0;

		math::array<double, 1> vector = {0.0, 1.0, };
		
		solvers::least_squares(matrix, vector);

		CHECK(vector[0] == 0.0_a);
		CHECK(vector[1] == 0.5_a);

  }

	SECTION("Symmetric real 2x2 -- 1"){
	
		using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 0.89653;
		matrix[0][1] = 0.41072;
		matrix[1][0] = 0.41072;
		matrix[1][1] = 0.69479;

		math::array<double, 1> vector = {0.21563, 0.40103, };
		
		solvers::least_squares(matrix, vector);

		CHECK(vector[0] == -0.0327899021_a);
		CHECK(vector[1] == 0.596579_a);
		
  }
	
	SECTION("Symmetric real 2x2 -- 2"){
	
		using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 6432.12;
		matrix[0][1] = 4502.48;
		matrix[1][0] = 4502.48;
		matrix[1][1] = 3151.74;

		math::array<double, 1> vector({1.0, 1.0});
		
		solvers::least_squares(matrix, vector);

		CHECK(vector[0] == -30.882245351_a);
		CHECK(vector[1] ==  44.1177546523_a);
		
  }

}

#endif

#endif
