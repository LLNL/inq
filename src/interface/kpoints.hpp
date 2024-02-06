/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__KPOINTS
#define INQ__INTERFACE__KPOINTS

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <ions/brillouin.hpp>

namespace inq {
namespace interface {

struct {

	std::string name() const {
		return "kpoints";
	}

	std::string one_line() const {
		return "Specifies the kpoints used to sample the Brillouin zone in the simulation";
	}
	
	void help() const {
		
		std::cout << R""""(

The 'kpoints' command
==================

These are the options available:

- `kpoints`

   Prints a list of the current kpoints in the simulation.

   Example: `inq kpoints`

- `kpoints gamma`

   Example: `inq kpoints gamma`


- `kpoints grid <nx> <ny> <nz> [shifted]`

   Examples: `inq kpoints grid 8 8 8`
             `inq kpoints grid 4 4 4 shifted`

)"""";
		
		exit(0);
	}

  void operator()() const {
		auto bz = ions::brillouin::load(".inq/default_brillouin");
		if(input::environment::global().comm().root()) std::cout << bz;
	}

  void gamma() const {
		auto bz = ions::brillouin(systems::ions::load(".inq/default_ions"), input::kpoints::gamma());
		bz.save(input::environment::global().comm(), ".inq/default_brillouin");
	}
  
	template <typename ArgsType>
	void command(ArgsType const & args, bool const quiet) const {

    if(args.size() == 0) {
			operator()();
			exit(0);
		}

    if(args.size() == 1 and args[0] == "gamma") {
			gamma();
      if(not quiet) operator()();
			exit(0);
		}
    
		std::cerr << "Error: Invalid syntax in the 'kpoints' command" << std::endl;
		exit(1);
	}
	
} const kpoints;

}
}
#endif

#ifdef INQ_INTERFACE_KPOINTS_UNIT_TEST
#undef INQ_INTERFACE_KPOINTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
