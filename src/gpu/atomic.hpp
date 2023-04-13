/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GPU__ATOMIC
#define INQ__GPU__ATOMIC

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/run.hpp>
#include <math/complex.hpp>
#include <math/vector3.hpp>

namespace inq {
namespace gpu {
namespace atomic {

template <typename Type1, typename Type2>
GPU_FUNCTION inline Type1 add(Type1 * val, Type2 const & incr){
#ifdef ENABLE_CUDA
#ifdef __CUDA_ARCH__
	static_assert(__CUDA_ARCH__ >= 600, "gpu library needs target gpu architecture >= 6.0");
#endif
  return atomicAdd(val, incr);
#else
  auto old_val = *val;
  *val += incr;
  return old_val;
#endif
}

GPU_FUNCTION inline complex add(complex * val, complex const & incr){
	auto re = add((double *) val, real(incr));
	auto im = add(((double *) val) + 1, imag(incr));
	return complex(re, im);
}

template <typename Type, typename Space>
GPU_FUNCTION inline auto add(vector3<Type, Space> * val, vector3<Type, Space> const & incr){
	auto v0 = add((Type *) val + 0, incr[0]) ;
	auto v1 = add((Type *) val + 1, incr[1]) ;
	auto v2 = add((Type *) val + 2, incr[2]) ;
	return vector3<Type, Space>{v0, v1, v2};
}

#ifdef ENABLE_CUDA
template <typename Type2>
GPU_FUNCTION inline long add(long * val, Type2 const & incr){
  return add((unsigned long long int *) val, incr);
}

template <typename Type2>
GPU_FUNCTION inline long add(size_t * val, Type2 const & incr){
	static_assert(sizeof(size_t) == sizeof(unsigned long long int));
  return add((unsigned long long int *) val, incr);
}
#endif

}
}
}
#endif

#ifdef INQ_GPU_ATOMIC_UNIT_TEST
#undef INQ_GPU_ATOMIC_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif

