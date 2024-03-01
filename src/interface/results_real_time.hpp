/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__RESULTS_REAL_TIME
#define INQ__INTERFACE__RESULTS_REAL_TIME

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <interface/actions.hpp>
#include <input/environment.hpp>
#include <real_time/results.hpp>

#include <utility>

namespace inq {
namespace interface {

struct {		

	std::string name() const {
		return "results real-time";
	}

	std::string one_line() const {
		return "Get information about the results obtained from a real-time calculation";
	}

	void help() const {
		
		std::cout << R""""(

The 'results real-time' command
==================

This command queries the results obtained from a real-time
calculation. Without arguments, it prints the values calculated.

The options allows you to query a specific value. In this case only
the value will be printed without any other text, so it is suitable
for easy parsing in scripting. The values are returned in atomic
units.

These are the available subcommands:

- `results real-time`

  When no arguments are given, print the values calculated.

  Example: `inq results real-time`.


- `results real-time total-steps`

  Returns the total number of real-time simulation steps done.

  Example: `inq results real-time total-steps`.


- `results real-time total-time`

  Returns the total simulated time (in atomic units).

  Example: `inq results real-time total-time`.

- `results real-time time [step]`

  Returns the time values. If not additional arguments are passed, inq
  prints the whole series for each time step. Alternatively, you can
  pass a step index to get the energy value.

  Note that for the moment inq uses a uniform time integration, so the
  time is just the step index times the time-step.

  Examples: `inq results real-time time`.
            `inq results real-time time 99`.


- `results real-time total-energy [step]`

  Returns the values of the total energy during the propagation. If
  not additional arguments are passed, inq prints the whole series for
  each time step in a two column format, the first column is the time
  and the second one is the energy. This output is suitable to view
  on a plotting program like gnuplot.

  Alternatively, you can pass a step index to get the energy value for
  that step.

  Examples: `inq results real-time total-energy`
            `inq results real-time total-energy 43`.


)"""";
	}

private:

	static auto load() {
		try { return real_time::results::load(".inq/default_results_real_time"); }
		catch(...){
			actions::error(input::environment::global().comm(), "Cannot find real-time results, run a real-time simulation first");
			exit(1);
		}
	}
	
public:
	
	void operator()() const {
		auto res = load();
		if(input::environment::global().comm().root()) std::cout << res;
	}

	auto total_steps() const {
		return load().total_steps;
	}
	
	auto total_time() const {
		return load().total_time;
	}

	auto time() const {
		return load().time;
	}
	
	auto total_energy() const {
		return load().total_energy;
	}

	auto dipole() const {
		auto && res = load();
		if(res.dipole.size() == 0)	actions::error(input::environment::global().comm(), "The dipole was not calculated during the real-time simulation");
		return res.dipole;
	}
	
	auto current() const {
		auto && res = load();
		if(res.current.size() == 0)	actions::error(input::environment::global().comm(), "The current was not calculated during the real-time simulation.");
		return res.current;
	}
	
private:

	template <typename ArgsType, typename ArrayType> 
	void array_output_vector(ArgsType args, ArrayType const & array, std::string const & command_name) const {
		
		if(args.size() == 1) {
			
			auto time_array = time();
			if(input::environment::global().comm().root()) {
				for(auto ii = 0ul; ii < time_array.size(); ii++) printf("%.20e\t%.20e\t%.20e\t%.20e\n", time_array[ii], array[ii][0], array[ii][1], array[ii][2]);
			}
			actions::normal_exit();
			
		} else if (args.size() == 2 or args.size() == 3) {
			auto index = utils::str_to<long>(args[1]);
			if(index < 0 or index >= (long) array.size()) actions::error(input::environment::global().comm(), "Invalid index ", index, " in the '" + command_name + "' command");
			
			if(args.size() == 2) {
				if(input::environment::global().comm().root()) printf("%.20e\t%.20e\t%.20e\n", array[index][0], array[index][1], array[index][2]);
				actions::normal_exit();
			}
			
			auto idir = utils::str_to_index(args[2]);
			if(idir == -1) actions::error(input::environment::global().comm(), "Invalid coordinate index in the '" + command_name + "' command");
			if(input::environment::global().comm().root()) printf("%.20e\n", array[index][idir]);
			actions::normal_exit();
			
		} else {
			actions::error(input::environment::global().comm(), "Invalid syntax in the '" + command_name + "' command");
		}
	}
	
public:
	
	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {

		if(args.size() == 0){
			operator()();
			actions::normal_exit();
		}

		if(args.size() == 1 and args[0] == "total-steps"){
			if(input::environment::global().comm().root()) printf("%ld\n", total_steps());
			actions::normal_exit();
		}
		
		if(args.size() == 1 and args[0] == "total-time"){
			if(input::environment::global().comm().root()) printf("%.20e\n", total_time());
			actions::normal_exit();
		}

		if(args[0] == "time"){
			auto time_array = time();
			if(args.size() == 1) {
				if(input::environment::global().comm().root()) {
					for(auto & ti : time_array) printf("%.20e\n", ti);
				}
			} else if (args.size() == 2) {
				if(input::environment::global().comm().root()) printf("%.20e\n", time_array[utils::str_to<long>(args[1])]);
			} else {
				actions::error(input::environment::global().comm(), "Invalid syntax in the 'results real-time time' command");
			}
			actions::normal_exit();
		}

		if(args[0] == "total-energy"){
			auto energy_array = total_energy();
			if(args.size() == 1) {

				auto time_array = time();
				if(input::environment::global().comm().root()) {
					for(auto ii = 0ul; ii < time_array.size(); ii++) printf("%.20e\t%.20e\n", time_array[ii], energy_array[ii]);
				}
				
			} else if (args.size() == 2) {
				if(input::environment::global().comm().root()) printf("%.20e\n", energy_array[utils::str_to<long>(args[1])]);
			} else {
				actions::error(input::environment::global().comm(), "Invalid syntax in the 'results real-time total-energy' command");
			}
			actions::normal_exit();
		}

		if(args[0] == "dipole") array_output_vector(args, dipole(), "result real-time dipole");
		if(args[0] == "current") array_output_vector(args, current(), "result real-time current");
				
		actions::error(input::environment::global().comm(), "Invalid syntax in the 'results real-time' command");
	}
	
} const results_real_time;

}
}
#endif

#ifdef INQ_INTERFACE_RESULTS_REAL_TIME_UNIT_TEST
#undef INQ_INTERFACE_RESULTS_REAL_TIME_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
