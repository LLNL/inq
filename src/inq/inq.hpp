/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INQ
#define INQ__INQ

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <bomd/propagate.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <operations/io.hpp>
#include <real_time/propagate.hpp>
#include <utils/match.hpp>
#include <perturbations/kick.hpp>
#include <perturbations/laser.hpp>
#include <perturbations/sum.hpp>
#include <observables/spectrum.hpp>
#endif

#ifdef INQ_INQ_INQ_UNIT_TEST
#undef INQ_INQ_INQ_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
