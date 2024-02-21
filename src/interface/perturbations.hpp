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

Note that the perturbations can be applied in different ways depending
on the periodicity of the system. The gauge selection is done
automatically by the code.

These are the uses for the command:

- `perturbations`

  Without any arguments, `perturbations` prints a list of the perturbations that are currently defined.

  Example: `inq perturbations`.


- `perturbations clear`

  Removes any perturbations present in the system.

  Example: `inq perturbations clear`


- `perturbations kick <px> <py> <pz>`

  Adds a kick perturbation. This corresponds of a uniform
  time-dependent electric field. The time profile is a delta function
  at time 0. The direction and intensity of the field is given by the
  values _px_, _py_ and _pz_.

  A kick has a flat profile in frequency space, so it will excite all
  frequencies equally. This is useful to obtain linear response
  properties.

  For finite systems this kick is applied in the length gauge as a
  phase to the orbitals (this avoid time discretization errors). For
  periodic systems, the kick is applied in the velocity gauge as a
  uniform vector potential.

  Example: `inq perturbations kick 0.0 0.0 0.1`


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

	void kick(vector3<double, cartesian> const & polarization) const {
		auto perturbations = perturbations::blend::load(".inq/default_perturbations");
		auto cell = systems::ions::load(".inq/default_ions").cell();
		perturbations.add(perturbations::kick(cell, polarization));
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

		if(args[0] == "kick"){

			if(args.size() != 4) {
				if(input::environment::global().comm().root()) std::cerr << "Error: Invalid arguments for 'perturbations kick'." << std::endl;
				exit(1);
			}

			kick({str_to<double>(args[1]), str_to<double>(args[2]), str_to<double>(args[3])});
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
