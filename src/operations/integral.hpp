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
#include <operations/sum.hpp>

namespace operations {

  template <class field_type>
  auto integral(const field_type & phi){
		auto integral_value = phi.basis().volume_element()*sum(phi.linear());

		if(phi.basis().dist().parallel()){
			phi.basis().dist().comm().all_reduce_in_place_n(&integral_value, 1, std::plus<>{});
		}

		return integral_value;
	}

  template <class field_type, class binary_op>
  auto integral(const field_type & phi1, const field_type & phi2, const binary_op op){
		assert(phi1.basis() == phi2.basis());

		auto integral_value = phi1.basis().volume_element()*operations::sum(phi1.linear(), phi2.linear(), op);
		
		if(phi1.basis().dist().parallel()){
			phi1.basis().dist().comm().all_reduce_in_place_n(&integral_value, 1, std::plus<>{});
		}

		return integral_value;
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

	auto comm = boost::mpi3::environment::get_world_instance();
		
	basis::trivial bas(N, comm);
	
	SECTION("Integral double"){
		
		basis::field<basis::trivial, double> aa(bas);

		aa = 1.0;

		REQUIRE(operations::integral(aa) == 1.0_a);

		for(int ii = 0; ii < aa.basis().dist().local_size(); ii++)	aa.linear()[ii] = aa.basis().dist().local_to_global(ii);

		REQUIRE(operations::integral(aa) == Approx(0.5*N*(N - 1.0)*bas.volume_element()));

	}

	SECTION("Integral complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		aa = complex(1.0, 1.0);

		REQUIRE(real(operations::integral(aa)) == 1.0_a);
		REQUIRE(imag(operations::integral(aa)) == 1.0_a);

		for(int ii = 0; ii < aa.basis().dist().local_size(); ii++) {
			auto iig = aa.basis().dist().local_to_global(ii);
			aa.linear()[ii] = complex(iig, -3.0*iig);
		}

		REQUIRE(real(operations::integral(aa)) == Approx(0.5*N*(N - 1.0)*bas.volume_element()));
		REQUIRE(imag(operations::integral(aa)) == Approx(-1.5*N*(N - 1.0)*bas.volume_element()));

	}

	SECTION("Integral product double"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		
		aa = 2.0;
		bb = 0.8;
		
		REQUIRE(operations::integral_product(aa, bb) == 1.6_a);
		
		for(int ii = 0; ii < aa.basis().dist().local_size(); ii++)	{
			auto iig = aa.basis().dist().local_to_global(ii);
			aa.linear()[ii] = pow(iig + 1, 2);
			bb.linear()[ii] = 1.0/(iig + 1);
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
		
		for(int ii = 0; ii < aa.basis().dist().local_size(); ii++)	{
			auto iig = aa.basis().dist().local_to_global(ii);
			aa.linear()[ii] = pow(iig + 1, 2)*exp(complex(0.0, M_PI/8 + M_PI/7*iig));
			bb.linear()[ii] = 1.0/(iig + 1)*exp(complex(0.0, M_PI/8 - M_PI/7*iig));
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
		for(int ii = 0; ii < aa.basis().dist().local_size(); ii++)	{
			auto iig = aa.basis().dist().local_to_global(ii);
			aa.linear()[ii] = sign*2.0*(iig + 1);
			bb.linear()[ii] = sign*1.0*(iig + 1);
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
		for(int ii = 0; ii < aa.basis().dist().local_size(); ii++)	{
			auto iig = aa.basis().dist().local_to_global(ii);
			aa.linear()[ii] = sign*2.0*(iig + 1)*exp(complex(0.0, 0.123*iig));
			bb.linear()[ii] = sign*1.0*(iig + 1)*exp(complex(0.0, 0.123*iig));
			sign *= -1.0;
		}
		
		REQUIRE(operations::integral_absdiff(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
		
}


#endif
#endif
