/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__REAL_TIME
#define INQ__INTERFACE__REAL_TIME

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
		return "real-time";
	}

	std::string one_line() const {
		return "Defines the parameters for a real-time calculation";
	}
	
	void help() const {
		
		std::cout << R""""(

The 'real-time' command
==================

This command defines the options for real-time self-consistency
calculations. These are the available options:

- `real-time`

  When no arguments are given, `real-time` will just print the
  currently defined options (including default values).

  Example: `inq real-time`.


- `real-time time-step <value> units`

  Sets the time step for the real-time integration. In most cases you
  want to pick the largest value that gives you a stable
  propagation. The default value is 0.01 atu.

  Example: 'inq real-time time-step 0.1 atu'.


- `real-time num-steps <value>`

  The number of time-steps in the time propagation. The default value
  is 100.

  Example: `inq real-time num-steps 10000`.


)"""";
	}

	void operator()() const {
		auto opts = options::real_time::load(".inq/default_real_time_options");
		if(input::environment::global().comm().root()) std::cout << opts;
	}

	void time_step(quantity<magnitude::time> dt) const {
		using namespace magnitude;

		auto opts = options::real_time::load(".inq/default_real_time_options").dt(dt);
		opts.save(input::environment::global().comm(), ".inq/default_real_time_options");
	}
	
	void num_steps(long nsteps) const {
		auto opts = options::real_time::load(".inq/default_real_time_options").num_steps(nsteps);
		opts.save(input::environment::global().comm(), ".inq/default_real_time_options");
	}

	void ions_static() const {
		auto opts = options::real_time::load(".inq/default_real_time_options").static_ions();
		opts.save(input::environment::global().comm(), ".inq/default_real_time_options");
	}

	void ions_impulsive() const {
		auto opts = options::real_time::load(".inq/default_real_time_options").impulsive();
		opts.save(input::environment::global().comm(), ".inq/default_real_time_options");
	}
	
	void ions_ehrenfest() const {
		auto opts = options::real_time::load(".inq/default_real_time_options").ehrenfest();
		opts.save(input::environment::global().comm(), ".inq/default_real_time_options");
	}
	
	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {
		using utils::str_to;
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}

		if(args.size() == 3 and (args[0] == "time-step")){
			time_step(str_to<double>(args[1])*magnitude::time::parse(args[2]));
			if(not quiet) operator()();
			exit(0);
		}
		
		if(args.size() == 2 and (args[0] == "num-steps")){
			num_steps(str_to<long>(args[1]));
			if(not quiet) operator()();
			exit(0);
		}

		if(args[0] == "ions"){
			args.erase(args.begin());

			if(args.size() == 1 and args[0] == "static"){
				ions_static();
				if(not quiet) operator()();
				exit(0);
			}
			
			if(args.size() == 1 and args[0] == "impulsive"){
				ions_impulsive();
				if(not quiet) operator()();
				exit(0);
			}
			
			if(args.size() == 1 and args[0] == "ehrenfest"){
				ions_ehrenfest();
				if(not quiet) operator()();
				exit(0);
			}
			
			if(input::environment::global().comm().root()) std::cerr << "Error: Invalid arguments for 'real-time ions' command" << std::endl;
			exit(1);
		}
				
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in 'real-time' command" << std::endl;
		exit(1);
	}
	
} const real_time;

}
}
#endif

#ifdef INQ_INTERFACE_REAL_TIME_UNIT_TEST
#undef INQ_INTERFACE_REAL_TIME_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
