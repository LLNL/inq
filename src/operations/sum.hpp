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

#include <gpu/reduce.hpp>
#include <math/array.hpp>
#include <math/zero.hpp>

#include <cassert>
#include <functional> //std::plus
#include <numeric>

#include <utils/profiling.hpp>

namespace inq {
namespace operations {

template <class array_type>
typename array_type::element sum(const array_type & phi){

	CALI_CXX_MARK_SCOPE("sum(1arg)");
	if(phi.size() == 0) return zero<typename array_type::element>();
	return gpu::run(gpu::reduce(phi.size()), gpu::array_access<decltype(begin(phi))>{begin(phi)});
}

template <class array1_type, class array2_type, class binary_op>
auto sum(const array1_type & phi1, const array2_type & phi2, const binary_op op){

	CALI_CXX_MARK_SCOPE("sum(2arg)");

	assert(phi1.size() == phi2.size());
	
	using return_type = decltype(op(phi1[0], phi2[0]));

	return_type initial = 0.0;
	return std::inner_product(phi1.begin(), phi1.end(), phi2.begin(), initial, std::plus<>(), op);
}

template <class ArrayType1, class ArrayType2>
auto sum_product(ArrayType1 & phi1, ArrayType2 const & phi2){
	return sum(phi1, phi2, std::multiplies<>());
}
	
}
}
#endif

#ifdef INQ_OPERATIONS_SUM_UNIT_TEST
#undef INQ_OPERATIONS_SUM_UNIT_TEST

#include <math/complex.hpp>
#include <basis/field.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("function operations::sum", "[operations::sum]") {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
	
	const int N = 21111;

	SECTION("Sum double"){

		math::array<double, 1> aa(N);

		aa.fill(1.0);

		CHECK(operations::sum(aa) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa[ii] = ii;

		CHECK(operations::sum(aa) == Approx(0.5*N*(N - 1.0)));

	}
	
	SECTION("Sum complex"){
		
		math::array<complex, 1> aa(N);
		
		aa.fill(complex(1.0, 1.0));

		CHECK(real(operations::sum(aa)) == Approx(N));
		CHECK(imag(operations::sum(aa)) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa[ii] = complex(ii, -3.0*ii);

		CHECK(real(operations::sum(aa)) == Approx(0.5*N*(N - 1.0)));
		CHECK(imag(operations::sum(aa)) == Approx(-1.5*N*(N - 1.0)));

	}

	SECTION("Sum product double"){

		math::array<double, 1> aa(N);
		math::array<double, 1> bb(N);

		aa.fill(2.0);
		bb.fill(0.8);
		
		CHECK(operations::sum_product(aa, bb) == Approx(1.6*N));
		
		for(int ii = 0; ii < N; ii++)	{
			aa[ii] = pow(ii + 1, 2);
			bb[ii] = 1.0/(ii + 1);
		}
		
		CHECK(operations::sum_product(aa, bb) == Approx(0.5*N*(N + 1.0)));
		
	}
	
	SECTION("Sum product complex"){

		math::array<complex, 1> aa(N);
		math::array<complex, 1> bb(N);
		
		aa.fill(complex(2.0, -0.3));
		bb.fill(complex(0.8, 0.01));
		
		CHECK(real(operations::sum_product(aa, bb)) == Approx(1.603*N));
		CHECK(imag(operations::sum_product(aa, bb)) == Approx(-0.22*N));
		
		for(int ii = 0; ii < N; ii++)	{
			aa[ii] = pow(ii + 1, 2)*exp(complex(0.0, M_PI/8 + M_PI/7*ii));
			bb[ii] = 1.0/(ii + 1)*exp(complex(0.0, M_PI/8 - M_PI/7*ii));
		}
		
		CHECK(real(operations::sum_product(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)));
		CHECK(real(operations::sum_product(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)));
		
	}
}
#endif
