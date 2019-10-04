/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__ADD
#define OPERATIONS__ADD

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

#include <cstdlib>

namespace operations {

	template <class field_type>
	auto add(const field_type & t1, const field_type & t2){
		assert(t1.basis() == t2.basis());
		
		field_type tadd(t1.basis());

		//DATAOPERATIONS
		for(long ii = 0; ii < t1.basis().size(); ii++) tadd[ii] = t1[ii] + t2[ii];
		
		return tadd;
	}

	template <class field_type>
	auto add(const field_type & t1, const field_type & t2, const field_type & t3){
		assert(t1.basis() == t2.basis());
		assert(t1.basis() == t3.basis());
		
		field_type tadd(t1.basis());
		
		//DATAOPERATIONS
		for(long ii = 0; ii < t1.basis().size(); ii++) tadd[ii] = t1[ii] + t2[ii] + t3[ii];
		
		return tadd;
	}

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::add", "[operations::add]") {

	using namespace Catch::literals;

	SECTION("Add 2 double arrays"){
		
		const int N = 100;

		basis::trivial bas(N);
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		basis::field<basis::trivial, double> cc(bas);

		aa = 1.0;
		bb = 2.5;

		cc = operations::add(aa, bb);
		
		for(int ii = 0; ii < N; ii++) REQUIRE(cc[ii] == 3.5_a);

	}
	
	SECTION("Add 2 complex arrays"){
		
		const int N = 100;

		basis::trivial bas(N);
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		basis::field<basis::trivial, complex> cc(bas);

		aa = complex(1.0, -20.2);
		bb = complex(2.5, 1.2);

		cc = operations::add(aa, bb);
		
		for(int ii = 0; ii < N; ii++){
			REQUIRE(real(cc[ii]) == 3.5_a);
			REQUIRE(imag(cc[ii]) == -19.0_a);
		}

	}

	SECTION("Add 3 double arrays"){
		
		const int N = 100;

		basis::trivial bas(N);
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		basis::field<basis::trivial, double> cc(bas);
		basis::field<basis::trivial, double> dd(bas);

		aa = 1.0;
		bb = 2.5;
		cc = -4.0;

		dd = operations::add(aa, bb, cc);
		
		for(int ii = 0; ii < N; ii++) REQUIRE(dd[ii] == -0.5_a);

	}
	
	SECTION("Add 3 complex arrays"){
		
		const int N = 100;

		basis::trivial bas(N);
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		basis::field<basis::trivial, complex> cc(bas);
		basis::field<basis::trivial, complex> dd(bas);

		aa = complex(1.0, -20.2);
		bb = complex(2.5, 1.2);
		cc = complex(-2.7, 8.6);

		dd = operations::add(aa, bb, cc);
		
		for(int ii = 0; ii < N; ii++){
			REQUIRE(real(dd[ii]) == 0.8_a);
			REQUIRE(imag(dd[ii]) == -10.4_a);
		}

	}	
}


#endif

#endif
