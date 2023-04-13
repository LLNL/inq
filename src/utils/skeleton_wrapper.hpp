/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__SKELETON_WRAPPER
#define INQ__UTILS__SKELETON_WRAPPER

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq {
namespace utils {

	template <class Base>
	struct skeleton_wrapper {
		skeleton_wrapper(Base const & bb)
			:base(bb){
		}

		Base const & base;
	};

}
}
#endif

#ifdef INQ_UTILS_SKELETON_WRAPPER_UNIT_TEST
#undef INQ_UTILS_SKELETON_WRAPPER_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif

