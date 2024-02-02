/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__GROUND_STATE
#define INQ__INTERFACE__GROUND_STATE

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <systems/electrons.hpp>

namespace inq {
namespace interface {

struct {
		
	std::string name() const {
		return "ground_state";
	}

	std::string one_line() const {
		return "Defines the how the ground_state is calculated.";
	}

	void operator()() const {
		auto gs_opts = options::ground_state::load(".inq/default_ground_state_options");
		std::cout << gs_opts;
	}

	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}
	}
	
} const ground_state;

}
}
#endif

#ifdef INQ_INTERFACE_GROUND_STATE_UNIT_TEST
#undef INQ_INTERFACE_GROUND_STATE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
