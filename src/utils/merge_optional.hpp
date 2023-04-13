/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__MERGE_OPTIONAL
#define INQ__UTILS__MERGE_OPTIONAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq {
namespace utils {

template <class opt_type>
opt_type merge_optional(const opt_type & option1, const opt_type & option2){
	if(option2) return option2;
	return option1;
}

}
}
#endif

#ifdef INQ_UTILS_MERGE_OPTIONAL_UNIT_TEST
#undef INQ_UTILS_MERGE_OPTIONAL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	CHECK(utils::merge_optional(std::optional<bool>{}, std::optional<bool>{}).has_value() == false);
	CHECK(utils::merge_optional(std::optional<int>{2}, std::optional<int>{}).value() == 2);
	CHECK(utils::merge_optional(std::optional<double>{}, std::optional<double>{3.4}).value() == 3.4_a);
	CHECK(utils::merge_optional(std::optional<std::string>{"uno"}, std::optional<std::string>{"dos"}).value() == "dos");
}
#endif
