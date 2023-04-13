/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PHYSICS__CONSTANTS
#define INQ__PHYSICS__CONSTANTS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq {
namespace physics {
namespace constants {

constexpr double proton_charge = -1.0;

}
}
}
#endif

#ifdef INQ_PHYSICS_CONSTANTS_UNIT_TEST
#undef INQ_PHYSICS_CONSTANTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif


