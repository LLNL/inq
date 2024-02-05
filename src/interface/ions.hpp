/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__IONS
#define INQ__INTERFACE__IONS

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <interface/cell.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace interface {

struct {
	
	std::string name() const {
		return "ions";
	}

	std::string one_line() const {
		return "Defines the ions in the simulation.";
	}

	void operator()() const {
		auto ions = systems::ions::load(".inq/default_ions");		
		if(input::environment::global().comm().root()) std::cout << ions;
	}

	void add(input::species const & sp, vector3<quantity<magnitude::length>> const & pos) const {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.insert(sp, pos);
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	void clear() const {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.clear();
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	void file(std::string const & filename) const {
		std::string extension = utils::lowercase(filename.substr(filename.find_last_of(".") + 1));

		if(extension == "xyz") {
			auto cell = systems::ions::load(".inq/default_ions").cell();
			auto ions = systems::ions::parse(filename, cell);
			ions.save(input::environment::global().comm(), ".inq/default_ions");
		} else {
			auto ions = systems::ions::parse(filename);
			ions.save(input::environment::global().comm(), ".inq/default_ions");
		}
	}
	
	void file(std::string const & filename, quantity<magnitude::length> const & radius) const {
		auto ions = systems::ions::parse(filename, radius);
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}

		if(args[0] == "clear"){

			if(args.size() != 1) {
				std::cerr << "Error: The 'ions clear' command doesn't take arguments." << std::endl;
				exit(1);
			}
			clear();
			if(not quiet) operator()();
			exit(0);
		}
 
		if(args[0] == "add"){

			if(args.size() != 6) {
				std::cerr << "Error: Wrong arguments for ions add.\nUse: inq ions add <symbol> <pos_x> <pos_y> <pos_z> <units>" << std::endl;
				exit(1);
			}

			auto symbol = args[1];
			auto units = magnitude::length::parse(args[5]);
			auto xx = atof(args[2].c_str())*units;
			auto yy = atof(args[3].c_str())*units;
			auto zz = atof(args[4].c_str())*units;
			
			add(symbol, {xx, yy, zz});
			if(not quiet) operator()();
			exit(0);
		}

		if(args.size() == 2 and args[0] == "file"){
			file(args[1]);
			if(not quiet) {
				interface::cell();
				operator()();
			}
			exit(0);
		}
		
		if(args.size() == 5 and args[0] == "file" and args[2] == "radius"){
			auto radius = atof(args[3].c_str())*magnitude::length::parse(args[4]);
			
			file(args[1], radius);
			if(not quiet) {
				interface::cell();
				operator()();
			}
			exit(0);
		}
		
		std::cerr << "Error: Invalid syntax in the ions command" << std::endl;
		exit(1);
	}
		
} const ions;

}
}
#endif

#ifdef INQ_INTERFACE_IONS_UNIT_TEST
#undef INQ_INTERFACE_IONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
