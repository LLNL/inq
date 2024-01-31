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

	void operator()() const {
		auto el_opts = options::electrons::load(".default_electrons_options");
		std::cout << el_opts;
	}

	void extra_states(int nstates) const{
		auto el_opts = options::electrons::load(".default_electrons_options").extra_states(nstates);
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}

	void extra_electrons(double nelectrons) const{
		auto el_opts = options::electrons::load(".default_electrons_options").extra_electrons(nelectrons);
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}
	
	void cutoff(quantity<magnitude::energy> ecut) const{
		auto el_opts = options::electrons::load(".default_electrons_options").cutoff(ecut);
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}

	void fourier_pseudo() const {
		auto el_opts = options::electrons::load(".default_electrons_options").fourier_pseudo();
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}

	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}
		
		if(args[0] == "extra_states"){

			if(args.size() == 1) {
				std::cerr << "Error: missing extra_states argument" << std::endl;
				exit(1);
			}

			if(args.size() >= 3) {
				std::cerr << "Error: too many arguments to extra_states argument" << std::endl;
				exit(1);
			}

			extra_states(atoi(args[1].c_str()));
			if(not quiet) operator()();
			exit(0);
		}
		
		if(args[0] == "extra_electrons"){

			if(args.size() == 1) {
				std::cerr << "Error: missing extra_electrons argument" << std::endl;
				exit(1);
			}

			if(args.size() >= 3) {
				std::cerr << "Error: too many arguments to extra_electrons argument" << std::endl;
				exit(1);
			}

			extra_electrons(atof(args[1].c_str()));
			if(not quiet) operator()();
			exit(0);
		}
		
		std::cerr << "Invalid syntax in 'electrons' command" << std::endl;
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
