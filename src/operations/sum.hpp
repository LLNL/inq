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
		//DATAOPERATIONS
		return std::accumulate(phi.begin(), phi.end(), (typename array_type::value_type) 0.0);
	}
	
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::sum", "[operations::sum]") {

	using namespace Catch::literals;
	
	const int N = 1000;
	
	basis::trivial bas(N);
	
	SECTION("Sum double"){
		
		basis::field<basis::trivial, double> aa(bas);

		aa = 1.0;

		REQUIRE(operations::sum(aa) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa[ii] = ii;

		REQUIRE(operations::sum(aa) == Approx(0.5*N*(N - 1.0)));

	}
	
	SECTION("Sum complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		aa = complex(1.0, 1.0);

		REQUIRE(real(operations::sum(aa)) == Approx(N));
		REQUIRE(imag(operations::sum(aa)) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa[ii] = complex(ii, -3.0*ii);

		REQUIRE(real(operations::sum(aa)) == Approx(0.5*N*(N - 1.0)));
		REQUIRE(imag(operations::sum(aa)) == Approx(-1.5*N*(N - 1.0)));

	}

}


#endif
#endif
