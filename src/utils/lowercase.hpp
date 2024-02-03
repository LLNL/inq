/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__LOWERCASE
#define INQ__UTILS__LOWERCASE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq {
namespace utils {

auto lowercase(std::string str) {
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

}
}
#endif

#ifdef INQ_UTILS_LOWERCASE_UNIT_TEST
#undef INQ_UTILS_LOWERCASE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	CHECK(utils::lowercase("AAA") == "aaa");
	CHECK(utils::lowercase("casa") == "casa");	
	CHECK(utils::lowercase("1935") == "1935");	
	
}
#endif
