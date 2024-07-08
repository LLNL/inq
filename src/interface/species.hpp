/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__SPECIES
#define INQ__INTERFACE__SPECIES

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
		return "species";
	}

	constexpr auto one_line() const {
		return "Define the species used in the simulation";
	}
	
	constexpr auto help() const {
		return R""""(

The 'species' command
==================

This command defines speciess the actual simulation.

These are the options available:

-  Shell:  `species`
   Python: `species.status()`

   Shows the species currently defined in the simulation.

   Shell example:  `inq species`
   Python example: `pinq.species.status()`

-  Shell:  `species pseudo-set`
   Python: `species.pseudo_set()`

   Returns the current pseudopotential set that will be used by
   default for all species in the simulation.

   Shell example:  `inq species pseudo-set`
   Python example: `pinq.species.pseudo_set()`

-  Shell:  `species pseudo-set <set name>`
   Python: `species.pseudo_set("set_name")`

   Sets the pseudopotential set that will be used by default for all
   species in the simulation.

   Shell example:  `inq species pseudo-set sg15`
   Python example: `pinq.species.pseudo-set("sg15")`

-  Shell:  `species list-sets`
   Python: `species.list_sets()`

   Prints a list of the available pseudopotentials.

   Shell example:  `inq species list-sets`
   Python example: `pinq.species.list_sets()`

)"""";
	}
	
	static void status() {
		if(input::environment::global().comm().root()) std::cout << systems::ions::load(".inq/default_ions").species_list();
	}
	
	void operator()() const {
		status();
	}

	static void pseudo_set() {
		if(input::environment::global().comm().root()) std::cout << systems::ions::load(".inq/default_ions").species_list().pseudopotentials() << std::endl;
	}

	static void pseudo_set(std::string set_name) {
		auto ions = systems::ions::load(".inq/default_ions");

		std::replace(set_name.begin(), set_name.end(), '-', '_');
			
		std::stringstream ss;
		ss << set_name;
		try { ss >> ions.species_list().pseudopotentials(); }
		catch(...) {
			actions::error(input::environment::global().comm(), "Unknown pseudo-set '" + set_name + "'."); 
		}
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	static void list_sets() {
		auto list = pseudo::set_id::list();
		std::cout << "Available pseudopotential sets:" << std::endl;
		for(auto const & item : list) {
			std::cout << "  " << item << std::endl;
		}
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		using utils::str_to;
		
		if(args.size() == 0) {
			status();
			actions::normal_exit();
		}

		if(args.size() == 1 and args[0] == "pseudo-set") {
			pseudo_set();
			actions::normal_exit();
		}

		if(args.size() == 2 and args[0] == "pseudo-set") {
			pseudo_set(args[1]);
      actions::normal_exit();
    }

		if(args.size() == 1 and args[0] == "list-sets") {
			list_sets();
			actions::normal_exit();
		}
		
	}
	
#ifdef INQ_PYTHON_INTERFACE
	template <class PythonModule>
	void python_interface(PythonModule & module) const {
		namespace py = pybind11;
		using namespace pybind11::literals;
 
		auto sub = module.def_submodule(name(), help());
		sub.def("status", &status);
		sub.def("pseudo_set", &pseudo_set);
		sub.def("list_sets", &list_sets);
		
	}
#endif

} const species;

}
}
#endif

#ifdef INQ_INTERFACE_SPECIES_UNIT_TEST
#undef INQ_INTERFACE_SPECIES_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
