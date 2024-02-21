/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__PERTURBATIONS
#define INQ__INTERFACE__PERTURBATIONS

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <perturbations/blend.hpp>

namespace inq {
namespace interface {

struct {
	
	std::string name() const {
		return "perturbations";
	}

	std::string one_line() const {
		return "Defines the perturbations in a real-time simulation";
	}
	void help() const {
		
		std::cout << R""""(

The 'perturbations' command
==================

This command defines the perturbations applied to the system in a real-time simulations.

These are the uses for the command:

- `perturbations`

  Without any arguments, `perturbations` prints a list of the perturbations that are currently defined.

  Example: `inq perturbations`.


- `perturbations clear`

  Removes any perturbations present in the system.

  Example: `inq perturbations clear`


)"""";
	}
	
	void operator()() const {
		auto perturbations = perturbations::blend::load(".inq/default_perturbations");
		if(input::environment::global().comm().root()) std::cout << perturbations;
	}

	void clear() const {
		auto perturbations = perturbations::blend::load(".inq/default_perturbations");
		perturbations.clear();
		perturbations.save(input::environment::global().comm(), ".inq/default_perturbations");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

		using utils::str_to;
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}

		if(args[0] == "clear"){

			if(args.size() != 1) {
				if(input::environment::global().comm().root()) std::cerr << "Error: The 'perturbations clear' command doesn't take arguments." << std::endl;
				exit(1);
			}
			clear();
			if(not quiet) operator()();
			exit(0);
		}
		
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in the perturbations command" << std::endl;
		exit(1);
	}
		
} const perturbations;

}
}
#endif

#ifdef INQ_INTERFACE_PERTURBATIONS_UNIT_TEST
#undef INQ_INTERFACE_PERTURBATIONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
