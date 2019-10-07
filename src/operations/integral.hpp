/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__INTEGRAL
#define OPERATIONS__INTEGRAL

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

  template <class field_type>
  auto integral(const field_type & phi){
		//OPTIMIZATION we should use std::reduce here, but it is not available in C++14
		//DATAOPERATIONS
		return phi.basis().volume_element()*std::accumulate(phi.begin(), phi.end(), (typename field_type::value_type) 0.0);
	}

  template <class field_type, class binary_op>
  auto integral(const field_type & phi1, const field_type & phi2, const binary_op op){
		assert(phi1.basis() == phi2.basis());
		
		const typename field_type::value_type initial = 0.0;
		//DATAOPERATIONS
		return  phi1.basis().volume_element()*std::inner_product(phi1.begin(), phi1.end(), phi2.begin(), initial, std::plus<>(), op);
	}
	
  template <class field_type>
  auto integral_product(const field_type & phi1, const field_type & phi2){
		return integral(phi1, phi2, std::multiplies<>());
	}
	
  template <class field_type>
  auto integral_absdiff(const field_type & phi1, const field_type & phi2){
		return real(integral(phi1, phi2, [](auto t1, auto t2){return fabs(t1 - t2);}));
	}
	
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::integral", "[operations::integral]") {

	using namespace Catch::literals;
	
	const int N = 1000;
	
	basis::trivial bas(N);
	
	SECTION("Integral double"){
		
		basis::field<basis::trivial, double> aa(bas);

		aa = 1.0;

		REQUIRE(operations::integral(aa) == 1.0_a);

		for(int ii = 0; ii < N; ii++)	aa[ii] = ii;

		REQUIRE(operations::integral(aa) == Approx(0.5*N*(N - 1.0)*bas.volume_element()));

	}
	
	SECTION("Integral complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		aa = complex(1.0, 1.0);

		REQUIRE(real(operations::integral(aa)) == 1.0_a);
		REQUIRE(imag(operations::integral(aa)) == 1.0_a);

		for(int ii = 0; ii < N; ii++)	aa[ii] = complex(ii, -3.0*ii);

		REQUIRE(real(operations::integral(aa)) == Approx(0.5*N*(N - 1.0)*bas.volume_element()));
		REQUIRE(imag(operations::integral(aa)) == Approx(-1.5*N*(N - 1.0)*bas.volume_element()));

	}

	SECTION("Integral product double"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		
		aa = 2.0;
		bb = 0.8;
		
		REQUIRE(operations::integral_product(aa, bb) == 1.6_a);
		
		for(int ii = 0; ii < N; ii++)	{
			aa[ii] = pow(ii + 1, 2);
			bb[ii] = 1.0/(ii + 1);
		}
		
		REQUIRE(operations::integral_product(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
	
	SECTION("Integral product complex"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		
		aa = complex(2.0, -0.3);
		bb = complex(0.8, 0.01);
		
		REQUIRE(real(operations::integral_product(aa, bb)) == 1.603_a);
		REQUIRE(imag(operations::integral_product(aa, bb)) == -0.22_a);
		
		for(int ii = 0; ii < N; ii++)	{
			aa[ii] = pow(ii + 1, 2)*exp(complex(0.0, M_PI/8 + M_PI/7*ii));
			bb[ii] = 1.0/(ii + 1)*exp(complex(0.0, M_PI/8 - M_PI/7*ii));
		}
		
		REQUIRE(real(operations::integral_product(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
		REQUIRE(real(operations::integral_product(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
		
	}
	
	
	SECTION("Integral absdiff double"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		
		aa = -13.23;
		bb = -13.23;
		
		REQUIRE(fabs(operations::integral_absdiff(aa, bb)) < 1e-14);

		double sign = 1.0;
		for(int ii = 0; ii < N; ii++)	{
			aa[ii] = sign*2.0*(ii + 1);
			bb[ii] = sign*1.0*(ii + 1);
			sign *= -1.0;
		}
		
		REQUIRE(operations::integral_absdiff(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
	
	SECTION("Integral absdiff complex"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		
		aa = -13.23*exp(complex(0.0, M_PI/3.63));
		bb = -13.23*exp(complex(0.0, M_PI/3.63));
		
		REQUIRE(fabs(operations::integral_absdiff(aa, bb)) < 1e-14);

		double sign = 1.0;
		for(int ii = 0; ii < N; ii++)	{
			aa[ii] = sign*2.0*(ii + 1)*exp(complex(0.0, 0.123*ii));
			bb[ii] = sign*1.0*(ii + 1)*exp(complex(0.0, 0.123*ii));
			sign *= -1.0;
		}
		
		REQUIRE(operations::integral_absdiff(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
		
}


#endif
#endif
