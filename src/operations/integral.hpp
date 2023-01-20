/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__INTEGRAL
#define INQ__OPERATIONS__INTEGRAL

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

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <operations/sum.hpp>
#include <math/complex.hpp>

namespace inq {
namespace operations {

template <class BasisType, class ElementType>
auto integral(basis::field<BasisType, ElementType> const & phi){
	CALI_CXX_MARK_FUNCTION;
	
	auto integral_value = phi.basis().volume_element()*sum(phi.linear());
	if(phi.basis().comm().size() > 1) {
		phi.basis().comm().all_reduce_in_place_n(&integral_value, 1, std::plus<>{});
	}
	return integral_value;
}

template <class BasisType, class ElementType>
auto integral_sum(basis::field_set<BasisType, ElementType> const & phi){
	CALI_CXX_MARK_FUNCTION;
	
	auto integral_value = phi.basis().volume_element()*sum(phi.matrix().flatted());
	if(phi.full_comm().size() > 1) {
		phi.full_comm().all_reduce_in_place_n(&integral_value, 1, std::plus<>{});
	}
	return integral_value;
}

template <class BasisType, class ElementType1, class ElementType2, class BinaryOp>
auto integral(basis::field<BasisType, ElementType1> const & phi1, basis::field<BasisType, ElementType2> const & phi2, BinaryOp const op){
	CALI_CXX_MARK_FUNCTION;
	
	assert(phi1.basis() == phi2.basis());

	auto integral_value = phi1.basis().volume_element()*operations::sum(phi1.linear(), phi2.linear(), op);
	if(phi1.basis().comm().size() > 1) {
		phi1.basis().comm().all_reduce_in_place_n(&integral_value, 1, std::plus<>{});
	}
	return integral_value;
}

template <class BasisType, class ElementType1, class ElementType2, class BinaryOp>
auto integral_sum(basis::field_set<BasisType, ElementType1> const & phi1, basis::field_set<BasisType, ElementType2> const & phi2, BinaryOp const op){
	CALI_CXX_MARK_FUNCTION;
	
	assert(phi1.basis() == phi2.basis());

	auto integral_value = phi1.basis().volume_element()*operations::sum(phi1.matrix().flatted(), phi2.matrix().flatted(), op);
	if(phi1.full_comm().size() > 1) {
		phi1.full_comm().all_reduce_in_place_n(&integral_value, 1, std::plus<>{});
	}
	return integral_value;
}

template <class BasisType, class ElementType>
auto integral_abs(basis::field<BasisType, ElementType> const & phi){
	return integral(phi, phi, [](auto t1, auto t2){return fabs(t1);});
}

template <class BasisType, class ElementType1, class ElementType2>
auto integral_product(basis::field<BasisType, ElementType1> const & phi1, basis::field<BasisType, ElementType2> const & phi2){
	return integral(phi1, phi2, std::multiplies<>());
}

template <class BasisType, class ElementType1, class ElementType2>
auto integral_product_sum(basis::field_set<BasisType, ElementType1> const & phi1, basis::field_set<BasisType, ElementType2> const & phi2){
	return integral_sum(phi1, phi2, std::multiplies<>());
}

template <class BasisType, class ElementType1, class ElementType2>
auto integral_absdiff(basis::field<BasisType, ElementType1> const & phi1, basis::field<BasisType, ElementType2> const & phi2){
	return real(integral(phi1, phi2, [](auto t1, auto t2){return fabs(t1 - t2);}));
}

template <class BasisType, class ElementType1, class ElementType2>
auto integral_sum_absdiff(basis::field_set<BasisType, ElementType1> const & phi1, basis::field_set<BasisType, ElementType2> const & phi2){
	return real(integral_sum(phi1, phi2, [](auto t1, auto t2){return fabs(t1 - t2);}));
}

}
}

#ifdef INQ_OPERATIONS_INTEGRAL_UNIT_TEST
#undef INQ_OPERATIONS_INTEGRAL_UNIT_TEST

#include <basis/field.hpp>
#include <catch2/catch_all.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::integral", "[operations::integral]") {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int N = 1000;

	auto comm = boost::mpi3::environment::get_world_instance();
		
	basis::trivial bas(N, comm);
	
	SECTION("Integral double"){
		
		basis::field<basis::trivial, double> aa(bas);

		aa.fill(1.0);

		CHECK(operations::integral(aa) == 1.0_a);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++)	aa.linear()[ii] = aa.basis().part().local_to_global(ii).value();

		CHECK(operations::integral(aa) == Approx(0.5*N*(N - 1.0)*bas.volume_element()));

	}

	SECTION("Integral complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		aa.fill(complex(1.0, 1.0));

		CHECK(real(operations::integral(aa)) == 1.0_a);
		CHECK(imag(operations::integral(aa)) == 1.0_a);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) {
			auto iig = aa.basis().part().local_to_global(ii);
			aa.linear()[ii] = complex(iig.value(), -3.0*iig.value());
		}

		CHECK(real(operations::integral(aa)) == Approx(0.5*N*(N - 1.0)*bas.volume_element()));
		CHECK(imag(operations::integral(aa)) == Approx(-1.5*N*(N - 1.0)*bas.volume_element()));

	}

	SECTION("integral_sum double"){

		int nvec = 6;
		
		basis::field_set<basis::trivial, double> aa(bas, nvec);

		aa.fill(1.0);

		CHECK(operations::integral_sum(aa) == Approx(nvec));

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) {
			for(int ist = 0; ist < nvec; ist++){
				aa.matrix()[ii][ist] = aa.basis().part().local_to_global(ii).value();
			}
		}

		CHECK(operations::integral_sum(aa) == Approx(nvec*0.5*N*(N - 1.0)*bas.volume_element()));
	}
	
	SECTION("Integral product double"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		
		aa.fill(2.0);
		bb.fill(0.8);
		
		CHECK(operations::integral_product(aa, bb) == 1.6_a);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++)	{
			auto iig = aa.basis().part().local_to_global(ii);
			aa.linear()[ii] = pow(iig.value() + 1, 2);
			bb.linear()[ii] = 1.0/(iig.value() + 1);
		}
		
		CHECK(operations::integral_product(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
	
	SECTION("Integral product complex"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		
		aa.fill(complex(2.0, -0.3));
		bb.fill(complex(0.8, 0.01));
		
		CHECK(real(operations::integral_product(aa, bb)) == 1.603_a);
		CHECK(imag(operations::integral_product(aa, bb)) == -0.22_a);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++)	{
			auto iig = aa.basis().part().local_to_global(ii);
			aa.linear()[ii] = pow(iig.value() + 1, 2)*exp(complex(0.0, M_PI/8 + M_PI/7*iig.value()));
			bb.linear()[ii] = 1.0/(iig.value() + 1)*exp(complex(0.0, M_PI/8 - M_PI/7*iig.value()));
		}
		
		CHECK(real(operations::integral_product(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
		CHECK(real(operations::integral_product(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
		
	}
	
	
	SECTION("Integral absdiff double"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		
		aa.fill(-13.23);
		bb.fill(-13.23);
		
		CHECK(fabs(operations::integral_absdiff(aa, bb)) < 1e-14);

		double sign = 1.0;
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++)	{
			auto iig = aa.basis().part().local_to_global(ii);
			aa.linear()[ii] = sign*2.0*(iig.value() + 1);
			bb.linear()[ii] = sign*1.0*(iig.value() + 1);
			sign *= -1.0;
		}
		
		CHECK(operations::integral_absdiff(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
	
	SECTION("Integral absdiff complex"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		
		aa.fill(-13.23*exp(complex(0.0, M_PI/3.63)));
		bb.fill(-13.23*exp(complex(0.0, M_PI/3.63)));
		
		CHECK(fabs(operations::integral_absdiff(aa, bb)) < 1e-14);

		double sign = 1.0;
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++)	{
			auto iig = aa.basis().part().local_to_global(ii);
			aa.linear()[ii] = sign*2.0*(iig.value() + 1)*exp(complex(0.0, 0.123*iig.value()));
			bb.linear()[ii] = sign*1.0*(iig.value() + 1)*exp(complex(0.0, 0.123*iig.value()));
			sign *= -1.0;
		}
		
		CHECK(operations::integral_absdiff(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
		
	}

	SECTION("integral_sum_absdiff double"){
		int nvec = 6;
		
		basis::field_set<basis::trivial, double> aa(bas, nvec);
		basis::field_set<basis::trivial, double> bb(bas, nvec);
		
		aa.fill(-13.23);
		bb.fill(-13.23);
		
		CHECK(fabs(operations::integral_sum_absdiff(aa, bb)) < 1e-14);

		double sign = 1.0;
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++)	{
			for(int ist = 0; ist < nvec; ist++){
				auto iig = aa.basis().part().local_to_global(ii);
				aa.matrix()[ii][ist] = sign*2.0*(iig.value() + 1);
				bb.matrix()[ii][ist] = sign*1.0*(iig.value() + 1);
				sign *= -1.0;
			}
		}
		
		CHECK(operations::integral_sum_absdiff(aa, bb) == Approx(nvec*0.5*N*(N + 1.0)*bas.volume_element()));
		
	}
	
}


#endif
#endif
