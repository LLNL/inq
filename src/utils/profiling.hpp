/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__PROFILING
#define INQ__UTILS__PROFILING

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <kalimotxo/cali.h>
#endif

#ifdef INQ_UTILS_PROFILING_UNIT_TEST
#undef INQ_UTILS_PROFILING_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
