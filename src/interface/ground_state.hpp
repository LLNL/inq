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
		return "ground-state";
	}

	std::string one_line() const {
		return "Defines the how the ground_state is calculated";
	}
	
	void help() const {
		
		std::cout << R""""(

The 'ground-state' command
==================

This command defines the options for ground-state self-consistency
calculations. These are the available options:

- `ground-state`

  When no arguments are given, `ground-state` will just print the
  currently defined options (including default values).

  Example: `inq ground-state`.


- `ground-state max-steps <value>`

  Sets the maximum number of self-consistency steps that will be
  done. The default value is 200.

  Example: 'inq ground-state max-steps 250'.


- `ground-state tolerance <value>`

  The tolerance used to consider that the self-consistency iteration
  is converged. The default value is 1e-6.

  Example: `inq ground-state tolerance 1e-9`.


- `ground-state mixing <value>`

  Set the mixing factor for the self-consistency. The default value is 0.3.

  Example: `inq ground-state mixing 0.1`.


)"""";
	}

	void operator()() const {
		auto gs_opts = options::ground_state::load(".inq/default_ground_state_options");
		if(input::environment::global().comm().root()) std::cout << gs_opts;
	}

	void max_steps(int nsteps) const {
		using namespace magnitude;

		auto gs_opts = options::ground_state::load(".inq/default_ground_state_options").max_steps(nsteps);
		gs_opts.save(input::environment::global().comm(), ".inq/default_ground_state_options");
	}
	
	void tolerance(double tol) const {
		using namespace magnitude;

		auto gs_opts = options::ground_state::load(".inq/default_ground_state_options").energy_tolerance(tol*1.0_Ha);
		gs_opts.save(input::environment::global().comm(), ".inq/default_ground_state_options");
	}
	
	void mixing(double factor) const {
		using namespace magnitude;

		auto gs_opts = options::ground_state::load(".inq/default_ground_state_options").mixing(factor);
		gs_opts.save(input::environment::global().comm(), ".inq/default_ground_state_options");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		using utils::str_to;
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}

		if(args.size() == 2 and (args[0] == "max-steps")){
			max_steps(str_to<int>(args[1]));
			if(not quiet) operator()();
			exit(0);
		}
		
		if(args.size() == 2 and (args[0] == "tolerance")){
			tolerance(str_to<double>(args[1]));
			if(not quiet) operator()();
			exit(0);
		}

		if(args.size() == 2 and (args[0] == "mixing")){
			mixing(str_to<double>(args[1]));
			if(not quiet) operator()();
			exit(0);
		}
		
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in 'ground-state' command" << std::endl;
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
