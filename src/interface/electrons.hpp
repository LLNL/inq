/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__ELECTRONS
#define INQ__INTERFACE__ELECTRONS

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
		return "electrons";
	}

	std::string one_line() const {
		return "Defines the electrons in the simulation and how they are represented.";
	}

	void help() const {
		
		std::cout << R""""(

The 'electrons' command
==================

This command defines the electrons that are in the system and how they
are going to be represented through several values that can be set by
the user.

- `electrons`

  When no arguments are given, `electrons` will just print the
  currently defined options (including default values).

  Example: `inq electrons`.


- `electrons cutoff <value> <units>`

  Sets the energy cutoff for the simulation grid. A higher cutoff
  implies a more precise, but more costly, simulation. The value must
  be followed by its units, check `inq help units` for details on what
  units are available.

  Example: `inq electrons cutoff 30.0 Ry`.


- `electrons spin <value>`

  Sets the spin configuration used in the simulation. The valid values
  are 'unpolarized' (the default), 'polarized' and 'non-collinear'.

  Example: `inq electrons spin polarized`.


- `electrons extra-electrons <value>`

  Inq determines the number of electrons from the ions present in the
  system. Using this variable you can add or remove electrons from the
  simulation. Positive values add electrons and negative values remove
  them. The number of electrons does not have to be an integer, you
  can add fractions of an electrons (given as a decimal number).

  Note that in first principles simulations the electrons are not
  associated a priory to a certain atom. So when you add an electron
  there is no concept of 'where' you put it. This will be determined
  by the ground-state optimization.

  Example: `inq electrons extra-electrons -0.5`.


- `electrons extra-states <value>`

  Inq automatically selects a number of states (orbitals, bands) that is
  enough to represent all the electrons in the system. In many cases
  you want additional states, and you do that using this command. The
  value must be a positive integer.

  Extra-states are necessary when setting an electronic temperature
  and to improve ground-state convergence.

  Example: `inq electrons extra-states 2`.


- `temperature <value> <units>`

  This command sets the temperature of the electrons in the
  ground-state optimization. The value must be positive and the units
  must be given. Check `inq help units` for details on what units are
  available. Most likely you want to use 'eV' or 'K'.

  Note that when you add a temperature you also need to specify
  extra-states.

  Example: `inq electrons temperature 273.15 Kelvin`.


)"""";
	}

	void operator()() const {
		auto el_opts = options::electrons::load(".inq/default_electrons_options");
		if(input::environment::global().comm().root()) std::cout << el_opts;
	}

	void extra_states(int nstates) const{
		auto el_opts = options::electrons::load(".inq/default_electrons_options").extra_states(nstates);
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}

	void extra_electrons(double nelectrons) const{
		auto el_opts = options::electrons::load(".inq/default_electrons_options").extra_electrons(nelectrons);
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}
	
	void cutoff(quantity<magnitude::energy> ecut) const{
		auto el_opts = options::electrons::load(".inq/default_electrons_options").cutoff(ecut);
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}

	void spacing(quantity<magnitude::length> const & val) const{
		auto el_opts = options::electrons::load(".inq/default_electrons_options").spacing(val);
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}
	
	void fourier_pseudo() const {
		auto el_opts = options::electrons::load(".inq/default_electrons_options").fourier_pseudo();
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}

	void spin_unpolarized() const {
		auto el_opts = options::electrons::load(".inq/default_electrons_options").spin_unpolarized();
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}

	void spin_polarized() const {
		auto el_opts = options::electrons::load(".inq/default_electrons_options").spin_polarized();
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}

	void spin_non_collinear() const {
		auto el_opts = options::electrons::load(".inq/default_electrons_options").spin_non_collinear();
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}

	void temperature(quantity<magnitude::energy> temp) const{
		auto el_opts = options::electrons::load(".inq/default_electrons_options").temperature(temp);
		el_opts.save(input::environment::global().comm(), ".inq/default_electrons_options");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

		using utils::str_to;
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}
		
		if(args[0] == "extra-states"){

			if(args.size() == 1) {
				if(input::environment::global().comm().root()) std::cerr << "Error: missing extra_states argument" << std::endl;
				exit(1);
			}

			if(args.size() >= 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: too many arguments to extra_states argument" << std::endl;
				exit(1);
			}

			extra_states(str_to<int>(args[1]));
			if(not quiet) operator()();
			exit(0);
		}
		
		if(args[0] == "extra-electrons"){

			if(args.size() == 1) {
				if(input::environment::global().comm().root()) std::cerr << "Error: missing extra_electrons argument" << std::endl;
				exit(1);
			}

			if(args.size() >= 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: too many arguments to extra_electrons argument" << std::endl;
				exit(1);
			}

			extra_electrons(str_to<double>(args[1]));
			if(not quiet) operator()();
			exit(0);
		}

		if(args[0] == "cutoff"){

			if(args.size() < 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: missing cutoff arguments. Use 'cutoff <value> <units>'" << std::endl;
				exit(1);
			}

			if(args.size() > 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: too many arguments to cutoff argument" << std::endl;
				exit(1);
			}

			cutoff(magnitude::energy::parse(str_to<double>(args[1]), args[2]));
			
			if(not quiet) operator()();
			exit(0);
		}

		if(args[0] == "spacing"){

			if(args.size() < 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: missing spacing arguments. Use 'spacing <value> <units>'" << std::endl;
				exit(1);
			}

			if(args.size() > 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: too many arguments to spacing argument" << std::endl;
				exit(1);
			}

			spacing(str_to<double>(args[1])*magnitude::length::parse(args[2]));
			
			if(not quiet) operator()();
			exit(0);
		}
		
		if(args.size() == 2 and args[0] == "spin" and args[1] == "unpolarized"){
			spin_unpolarized();
			if(not quiet) operator()();
			exit(0);
		}

		if(args.size() == 2 and args[0] == "spin" and args[1] == "polarized"){
			spin_polarized();
			if(not quiet) operator()();
			exit(0);
		}

		if(args.size() == 2 and args[0] == "spin" and args[1] == "non-collinear") {
			spin_non_collinear();
			if(not quiet) operator()();
			exit(0);
		}

		if(args[0] == "temperature"){

			if(args.size() < 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: missing temperature arguments. Use 'temperature <value> <units>'." << std::endl;
				exit(1);
			}

			if(args.size() > 3) {
				if(input::environment::global().comm().root()) std::cerr << "Error: too many arguments to temperature argument.  Use 'temperature <value> <units>'." << std::endl;
				exit(1);
			}

			temperature(magnitude::energy::parse(str_to<double>(args[1]), args[2]));
			
			if(not quiet) operator()();
			exit(0);
		}

		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in the 'electrons' command" << std::endl;
		exit(1);
	}
	
} const electrons;

}
}
#endif

#ifdef INQ_INTERFACE_ELECTRONS_UNIT_TEST
#undef INQ_INTERFACE_ELECTRONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
