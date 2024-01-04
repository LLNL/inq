/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__COMPLEX
#define INQ__MATH__COMPLEX

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <complex>
#include <gpu/run.hpp>

#ifdef ENABLE_GPU
#include <thrust/complex.h>
#endif

namespace inq {

#ifdef ENABLE_GPU
using complex = thrust::complex<double>;
using thrust::polar;
#else
using complex = std::complex<double>;
using std::polar;
#endif

GPU_FUNCTION inline double real(const double & x){
	return x;
}

GPU_FUNCTION inline double real(const complex & z){
	return z.real();
}

GPU_FUNCTION inline double imag(const double &){
	return 0.0;
}

GPU_FUNCTION inline auto imag(const complex & z){
	return z.imag();
}

GPU_FUNCTION inline auto norm(const double & x){
	return x*x;
}

GPU_FUNCTION inline double conj(const double & x){
	return x;
}

GPU_FUNCTION inline auto fabs(complex const & z){
	return abs(z);
}

}
#endif

#ifdef INQ_MATH_COMPLEX_UNIT_TEST
#undef INQ_MATH_COMPLEX_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;

	double xx = 203.42;

	CHECK(real(xx) == 203.42_a);
	CHECK(imag(xx) == 0.0_a);
	CHECK(norm(xx) == 41379.696_a);
	CHECK(conj(xx) == xx);

	complex zz{-654.21, 890.74};

	CHECK(real(zz) == -654.21_a);
	CHECK(imag(zz) == 890.74);
	CHECK(norm(zz) == 1221408.5_a);
	CHECK(fabs(zz) == 1105.1735_a);

}
#endif
