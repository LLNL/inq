/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__RAW_POINTER_CAST
#define INQ__UTILS__RAW_POINTER_CAST

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq{

template <typename Type>
auto raw_pointer_cast(Type * ptr){
	return ptr;
}

}
#endif

#ifdef INQ_UTILS_RAW_POINTER_CAST_UNIT_TEST
#undef INQ_UTILS_RAW_POINTER_CAST_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <math/complex.hpp>

using namespace inq;

TEMPLATE_TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG, int, double, int const, double const, complex, complex const) {

	TestType a = 1;
	
	CHECK(typeid(raw_pointer_cast(&a)) == typeid(TestType *));
	CHECK(raw_pointer_cast(&a) == &a);
	
}
#endif

