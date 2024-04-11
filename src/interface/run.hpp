/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__RUN
#define INQ__INTERFACE__RUN

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>
#include <real_time/results.hpp>

namespace inq {
namespace interface {

struct {

	constexpr auto name() const {
		return "run";
	}

	constexpr auto one_line() const {
		return "Runs the simulation";
	}
	
	constexpr auto help() const {
		
		return R""""(

The 'run' command
==================

This command runs the actual simulation. It requires all the
simulation parameters to be set before running.

Note that this is the expensive part in the execution on
inq. Depending on the system you are using you might want to execute
this in parallel or through a queuing system.

These are the options available:

- `run ground-state`

   Runs a ground-state calculation with fixed ions.

   Example: `inq run ground-state`


- `run real-time`

   Runs a real-time simulation.

   Example: `inq run real-time`


)"""";
	}

	void ground_state() const{
		auto ions = systems::ions::load(".inq/default_ions");

		auto bz = ions::brillouin(systems::ions::load(".inq/default_ions"), input::kpoints::gamma());

		try { bz = ions::brillouin::load(".inq/default_brillouin"); }
		catch(...) {
			bz.save(input::environment::global().comm(), ".inq/default_brillouin");
		}
		
		systems::electrons electrons(ions, options::electrons::load(".inq/default_electrons_options"), bz);
 
		if(not electrons.try_load(".inq/default_orbitals")){
			ground_state::initial_guess(ions, electrons);
		}
		auto result = ground_state::calculate(ions, electrons, options::theory::load(".inq/default_theory"), options::ground_state::load(".inq/default_ground_state_options").calculate_forces());

		result.save(input::environment::global().comm(), ".inq/default_results_ground_state");
		electrons.save(".inq/default_orbitals");
	}

	void real_time() const{
		auto ions = systems::ions::load(".inq/default_ions");

		auto bz = ions::brillouin(systems::ions::load(".inq/default_ions"), input::kpoints::gamma());

		try { bz = ions::brillouin::load(".inq/default_brillouin"); }
		catch(...) {
			bz.save(input::environment::global().comm(), ".inq/default_brillouin");
		}
		
		systems::electrons electrons(ions, options::electrons::load(".inq/default_electrons_options"), bz);
 
		if(not electrons.try_load(".inq/default_orbitals")) actions::error(input::environment::global().comm(), "Cannot load a ground-state electron configuration for a real-time run.\n Please run a ground-state first.");

		auto opts = options::real_time::load(".inq/default_real_time_options");	
		auto res = real_time::results(".inq/default_results_real_time");
		res.obs = opts.observables_container();
		
		real_time::propagate(ions, electrons, [&res](auto obs){ res(obs); }, options::theory::load(".inq/default_theory"), opts, perturbations::blend::load(".inq/default_perturbations"));
		res.save(input::environment::global().comm());

	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		if(args.size() == 0) actions::error(input::environment::global().comm(), "Missing argument to the 'run' command");
		
		if(args.size() == 1 and args[0] == "ground-state") {
			ground_state();
			actions::normal_exit();
		}
				
		if(args.size() == 1 and args[0] == "real-time") {
			real_time();
			actions::normal_exit();
		}
		
		actions::error(input::environment::global().comm(), "Invalid syntax in the 'run' command");
	}
	
} const run;

}
}
#endif

#ifdef INQ_INTERFACE_RUN_UNIT_TEST
#undef INQ_INTERFACE_RUN_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
