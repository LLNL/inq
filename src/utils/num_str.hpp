/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__NUM_STR
#define INQ__UTILS__NUM_STR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq {
namespace util {

auto num_to_str(long num){
	char numcstr[12]; 
	snprintf(numcstr, 11, "%010ld", num);
	return std::string(numcstr);				
}


}
}
#endif

#ifdef INQ_UTILS_NUM_STR_UNIT_TEST
#undef INQ_UTILS_NUM_STR_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <math/complex.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}

#endif

