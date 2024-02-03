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

	void tolerance(double tol) const {
		using namespace magnitude;

		auto gs_opts = options::ground_state::load(".inq/default_ground_state_options").energy_tolerance(tol*1.0_Ha);
		gs_opts.save(input::environment::global().comm(), ".inq/default_ground_state_options");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}

		if(args.size() == 2 and (args[0] == "tolerance" or args[0] == "tol")){
			tolerance(atof(args[1].c_str()));
			if(not quiet) operator()();
			exit(0);
		}

		std::cerr << "Error: Invalid syntax in 'ground_state' command" << std::endl;
		exit(1);
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
