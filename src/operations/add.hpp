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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <gpu/run.hpp>
#include <algorithm>
#ifdef HAVE_CUDA
#include <thrust/transform.h>
#endif

namespace operations {

	/*

		Returns a field that has the sum of the values of t1 and t2.

	*/
	template <class field_type>
	auto add(const field_type & t1, const field_type & t2){
		assert(t1.basis() == t2.basis());
		
		field_type tadd(t1.basis());

		using type = typename field_type::element_type;

		//DATAOPERATIONS STL + THRUST TRANSFORM
#ifdef HAVE_CUDA
		thrust::transform(t1.linear().begin(), t1.linear().end(), t2.linear().begin(), tadd.linear().begin(), thrust::plus<type>());
#else
		std::transform(t1.linear().begin(), t1.linear().end(), t2.linear().begin(), tadd.linear().begin(), std::plus<type>());
#endif
		
		return tadd;
	}
	
	/*

		Returns a field that has the sum of the values of t1, t2 and t3.

	*/
	template <class field_type>
	field_type add(const field_type & t1, const field_type & t2, const field_type & t3){
		assert(t1.basis() == t2.basis());
		assert(t1.basis() == t3.basis());
		
		field_type tadd(t1.basis());

		//DATAOPERATIONS LOOP + GPU::RUN 1D
#ifdef HAVE_CUDA

		using type = typename field_type::element_type;

		auto t1p = t1.linear().begin();
		auto t2p = t2.linear().begin();
		auto t3p = t3.linear().begin();
		auto taddp = tadd.linear().begin();
		
		gpu::run(t1.basis().size(),
						 [=] __device__ (long ii){
							 taddp[ii] = t1p[ii] + t2p[ii] + t3p[ii];
						 });
#else
		for(long ii = 0; ii < t1.basis().size(); ii++) tadd.linear()[ii] = t1.linear()[ii] + t2.linear()[ii] + t3.linear()[ii];
#endif
		
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
		
		for(int ii = 0; ii < N; ii++) REQUIRE(cc.linear()[ii] == 3.5_a);

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
			REQUIRE(real(cc.linear()[ii]) == 3.5_a);
			REQUIRE(imag(cc.linear()[ii]) == -19.0_a);
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
		
		for(int ii = 0; ii < N; ii++) REQUIRE(dd.linear()[ii] == -0.5_a);

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
			REQUIRE(real(dd.linear()[ii]) == 0.8_a);
			REQUIRE(imag(dd.linear()[ii]) == -10.4_a);
		}

	}	
}


#endif

#endif
