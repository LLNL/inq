/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SUM
#define OPERATIONS__SUM

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


#include <cassert>
#include <numeric>

namespace operations {

  template <class array_type>
  auto sum(const array_type & phi){
		//OPTIMIZATION we should use std::reduce here, but it is not available in C++14
		//DATAOPERATIONS STL ACCUMULATE
		return std::accumulate(phi.begin(), phi.end(), (typename array_type::element_type) 0.0);
	}

  template <class array1_type, class array2_type, class binary_op>
  auto sum(const array1_type & phi1, const array2_type & phi2, const binary_op op){

		const typename array1_type::element_type initial = 0.0;
		//OPTIMIZATION we should use std::transform_reduce here, but it is not available in C++14
		//DATAOPERATIONS STL INNER_PRODUCT
		return std::inner_product(phi1.begin(), phi1.end(), phi2.begin(), initial, std::plus<>(), op);
	}
	
	template <class field_type>
  auto sum_product(const field_type & phi1, const field_type & phi2){
		return sum(phi1, phi2, std::multiplies<>());
	}
	
}

#ifdef UNIT_TEST

#include <basis/field.hpp>

#include <catch2/catch.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::sum", "[operations::sum]") {

	using namespace Catch::literals;
	
	const int N = 1111;
	
	basis::trivial bas(N);
	
	SECTION("Sum double"){
		
		basis::field<basis::trivial, double> aa(bas);

		aa = 1.0;

		REQUIRE(operations::sum(aa.linear()) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa.linear()[ii] = ii;

		REQUIRE(operations::sum(aa.linear()) == Approx(0.5*N*(N - 1.0)));

	}
	
	SECTION("Sum complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		aa = complex(1.0, 1.0);

		REQUIRE(real(operations::sum(aa.linear())) == Approx(N));
		REQUIRE(imag(operations::sum(aa.linear())) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa.linear()[ii] = complex(ii, -3.0*ii);

		REQUIRE(real(operations::sum(aa.linear())) == Approx(0.5*N*(N - 1.0)));
		REQUIRE(imag(operations::sum(aa.linear())) == Approx(-1.5*N*(N - 1.0)));

	}

	SECTION("Sum product double"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		
		aa = 2.0;
		bb = 0.8;
		
		REQUIRE(operations::sum_product(aa.linear(), bb.linear()) == Approx(1.6*N));
		
		for(int ii = 0; ii < N; ii++)	{
			aa.linear()[ii] = pow(ii + 1, 2);
			bb.linear()[ii] = 1.0/(ii + 1);
		}
		
		REQUIRE(operations::sum_product(aa.linear(), bb.linear()) == Approx(0.5*N*(N + 1.0)));
		
	}
	
	SECTION("Sum product complex"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		
		aa = complex(2.0, -0.3);
		bb = complex(0.8, 0.01);
		
		REQUIRE(real(operations::sum_product(aa.linear(), bb.linear())) == Approx(1.603*N));
		REQUIRE(imag(operations::sum_product(aa.linear(), bb.linear())) == Approx(-0.22*N));
		
		for(int ii = 0; ii < N; ii++)	{
			aa.linear()[ii] = pow(ii + 1, 2)*exp(complex(0.0, M_PI/8 + M_PI/7*ii));
			bb.linear()[ii] = 1.0/(ii + 1)*exp(complex(0.0, M_PI/8 - M_PI/7*ii));
		}
		
		REQUIRE(real(operations::sum_product(aa.linear(), bb.linear())) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)));
		REQUIRE(real(operations::sum_product(aa.linear(), bb.linear())) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)));
		
	}
	
	
}


#endif
#endif
