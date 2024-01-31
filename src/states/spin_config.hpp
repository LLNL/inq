/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__STATES__SPIN_CONFIG
#define INQ__STATES__SPIN_CONFIG

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/complex.hpp>
#include <gpu/array.hpp>
#include <operations/sum.hpp>
#include <parallel/partition.hpp>

namespace inq {
namespace states {

enum class spin_config { UNPOLARIZED, POLARIZED, NON_COLLINEAR };

}
}
#endif

#ifdef INQ_STATES_SPIN_CONFIG_UNIT_TEST
#undef INQ_STATES_SPIN_CONFIG_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace Catch::literals;
	using namespace inq;

}
#endif
