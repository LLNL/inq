/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__LEAST_SQUARES
#define INQ__SOLVERS__LEAST_SQUARES

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


#include <math/array.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <FC.h>

#include <stdexcept>
#include <tuple> //for get

#define dgelss FC_GLOBAL(dgelss, DGELSS) 
extern "C" void dgelss(const int & m, const int & n, const int & nrhs, double * a, const int & lda, double * b, const int & ldb, double * s, const double & rcond, int & rank, double * work, const int & lwork, int & info);

namespace inq {
namespace solvers {

template <class matrix_type, class vector_type>
void least_squares(matrix_type && matrix, vector_type & rhs){

	int mm = std::get<0>(sizes(matrix));
	int nn = std::get<1>(sizes(matrix));

	math::array<double, 1> ss(mm);

	int rank, info;
	double dwork;
		
	//DATAOPERATIONS RAWLAPACK dgelss
	dgelss(mm, nn, 1, raw_pointer_cast(matrix.data_elements()), mm, raw_pointer_cast(rhs.data_elements()), mm, raw_pointer_cast(ss.base()), -1.0, rank, &dwork, -1, info);

	if(info != 0) std::runtime_error("inq error: dgelss in least_squares failed with info = " + std::to_string(info));
	
	auto work = (double *) malloc(int(dwork)*sizeof(double));

	dgelss(mm, nn, 1, raw_pointer_cast(matrix.data_elements()), mm, raw_pointer_cast(rhs.data_elements()), mm, raw_pointer_cast(ss.data_elements()), -1.0, rank, work, int(dwork), info);

	if(info != 0) std::runtime_error("inq error: dgelss in least_squares failed with info = " + std::to_string(info));
		
	free(work);
		
}

}
}

///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_LEAST_SQUARES_UNIT_TEST
#undef INQ_SOLVERS_LEAST_SQUARES_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("function solvers::least_squares", "[solvers::least_squares]") {

	SECTION("Diagonal real 2x2"){
	
		using namespace inq;
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
	
		using namespace inq;
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
	
		using namespace inq;
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
