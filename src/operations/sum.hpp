/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SUM
#define OPERATIONS__SUM

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/reduce.hpp>
#include <gpu/array.hpp>
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
	auto init = zero<typename array_type::element>();
	if(phi.size() == 0) return init;
	return gpu::run(gpu::reduce(phi.size()), init, [ph = begin(phi)] GPU_LAMBDA (auto ii) { return ph[ii]; });
}

template <class ArrayType, class UnaryOp>
auto sum(ArrayType const & arr, UnaryOp const & op){

	CALI_CXX_MARK_SCOPE("sum(1arg, op)");

	using return_type = decltype(op(arr[0]));

	return_type val = 0.0;
	for(int ii = 0; ii < arr.size(); ii++) val += op(arr[ii]);
	return val;
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
auto sum_product(ArrayType1 const & phi1, ArrayType2 const & phi2){
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

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
	
	const int N = 21111;

	SECTION("Sum double"){

		gpu::array<double, 1> aa(N);

		aa.fill(1.0);

		CHECK(operations::sum(aa) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa[ii] = ii;

		CHECK(operations::sum(aa) == Approx(0.5*N*(N - 1.0)));
		CHECK(operations::sum(aa, [](auto xx){ return -2.0*xx; }) == Approx(-N*(N - 1.0)));
	}
	
	SECTION("Sum complex"){
		
		gpu::array<complex, 1> aa(N);
		
		aa.fill(complex(1.0, 1.0));

		CHECK(real(operations::sum(aa)) == Approx(N));
		CHECK(imag(operations::sum(aa)) == Approx(N));

		for(int ii = 0; ii < N; ii++)	aa[ii] = complex(ii, -3.0*ii);

		CHECK(real(operations::sum(aa)) == Approx(0.5*N*(N - 1.0)));
		CHECK(imag(operations::sum(aa)) == Approx(-1.5*N*(N - 1.0)));

		CHECK(real(operations::sum(aa, [](auto xx) { return conj(complex{0.0, -1.0}*xx); })) == Approx(-1.5*N*(N - 1.0)));
		CHECK(imag(operations::sum(aa, [](auto xx) { return conj(complex{0.0, -1.0}*xx); })) == Approx( 0.5*N*(N - 1.0)));
		
	}

	SECTION("Sum product double"){

		gpu::array<double, 1> aa(N);
		gpu::array<double, 1> bb(N);

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

		gpu::array<complex, 1> aa(N);
		gpu::array<complex, 1> bb(N);
		
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
