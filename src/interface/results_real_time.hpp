/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__RESULTS_REAL_TIME
#define INQ__INTERFACE__RESULTS_REAL_TIME

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <real_time/results.hpp>

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


)"""";
	}

private:

	static auto load() {
		try { return real_time::results::load(".inq/default_results_real_time"); }
		catch(...){
			if(input::environment::global().comm().root()) std::cerr << "Error: cannot find real-time results, run a real-time simulation first" << std::endl;
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
	
	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {

		if(args.size() == 0){
			operator()();
			exit(0);
		}

		if(args.size() == 1 and args[0] == "total-steps"){
			if(input::environment::global().comm().root()) printf("%ld\n", total_steps());
			exit(0);
		}
		
		if(args.size() == 1 and args[0] == "total-time"){
			if(input::environment::global().comm().root()) printf("%.20e\n", total_time());
			exit(0);
		}
					
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in the 'results real-time' command" << std::endl;
		exit(1);
    
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
