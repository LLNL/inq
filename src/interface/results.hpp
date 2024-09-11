/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__RESULTS
#define INQ__INTERFACE__RESULTS

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace interface {

struct {

	constexpr auto name() const {
		return "results";
	}

	constexpr auto one_line() const {
		return "Prints information about the available input results in inq";
	}
	
	constexpr auto help() const {
		
		return R""""(

Results
=======

The result command queries the results of a simulation that has been
run in inq. Specific versions of the command exist for each simulation
run. You can see the specific help for each one with:

Shell:
- `inq help results ground-state`
- `inq help results real-time`

Python:
- `help(pinq.results.ground_state)`
- `help(pinq.results.real_time)`

)"""";
	}
		
} const results ;

}
}
#endif

#ifdef INQ_INTERFACE_RESULTS_UNIT_TEST
#undef INQ_INTERFACE_RESULTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
