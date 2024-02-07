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
	void help() const {
		
		std::cout << R""""(

The 'ions' command
==================

This command defines the ions present in the system and its
coordinates. These can be given one by one using the `ions insert`
command or from a file using the `ions file` command.

If you are adding ions one by one or reading an XYZ file, you must
first declare the cell for system (see `inq help cell`).

These are the uses for the command:

- `ions`

  Without any arguments, `ions` prints a list of the ions currently in the system.

  Example: `inq ions`.


- `ions clear`

  Removes any ions present in the system. Does not change the cell.

  Example: `inq ions clear`


- `ions insert <symbol> <x> <y> <z> <units>`

  Add an ion of type _symbol_ at coordinates _x_, _y_, and _z_. The
  units for the coordinates must be specified (check `inq help units`
  for details).

  Example: `ions insert He  0.0 0. 0.0 2.0 angstrom'


- `ions insert fractional <symbol> <x> <y> <z>`

  Insert an ion of type _symbol_ at fractional coordinates _x_, _y_,
  and _z_.

  Example: `ions insert fractional Si 0.25 0.25 0.25'


- `ions file <file>`

  Read a coordinate file. The supported formats are POSCAR, CIF and
  XYZ. For POSCAR and CIF both the ions and the cell information will
  be read. For XYZ only the atomic positions are read, so the cell
  must be defined before reading the file.

  Example: 'inq ions file diamond.cif'.


- `ions file <file> radius <r> <units>`

  Reads a coordinate file and define a cell around the ions. The cell
  is orthorhombic and finite. The cell has the smallest possible size
  so that walls are at least distance 'r' from any atom. This is
  useful for molecular systems where you need to converge the size of
  the cell.  The units for the radius value must be specified (check
  `inq help units` for details).

  Note: right now this functionality is only implemented for XYZ files.

  For example 'inq ions file glucose.xyz radius 2.0 A'.


)"""";
	}
	
	void operator()() const {
		auto ions = systems::ions::load(".inq/default_ions");		
		if(input::environment::global().comm().root()) std::cout << ions;
	}

	void insert(input::species const & sp, vector3<quantity<magnitude::length>> const & pos) const {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.insert(sp, pos);
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	void insert_fractional(input::species const & sp, vector3<double, contravariant> const & pos) const {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.insert_fractional(sp, pos);
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
				if(input::environment::global().comm().root()) std::cerr << "Error: The 'ions clear' command doesn't take arguments." << std::endl;
				exit(1);
			}
			clear();
			if(not quiet) operator()();
			exit(0);
		}

 		if(args.size() >= 2 and args[0] == "insert" and args[1] == "fractional"){

			if(args.size() != 6) {
				if(input::environment::global().comm().root()) std::cerr << "Error: Wrong arguments for ions insert.\nUse: inq ions insert fractional <symbol> <x> <y> <z>" << std::endl;
				exit(1);
			}

			auto symbol = args[2];
			auto xx = atof(args[3].c_str());
			auto yy = atof(args[4].c_str());
			auto zz = atof(args[5].c_str());
			
			insert_fractional(symbol, {xx, yy, zz});
			if(not quiet) operator()();
			exit(0);
		}

		if(args[0] == "insert"){

			if(args.size() != 6) {
				if(input::environment::global().comm().root()) std::cerr << "Error: Wrong arguments for ions insert.\nUse: inq ions insert <symbol> <x> <y> <z> <units>" << std::endl;
				exit(1);
			}

			auto symbol = args[1];
			auto units = magnitude::length::parse(args[5]);
			auto xx = atof(args[2].c_str())*units;
			auto yy = atof(args[3].c_str())*units;
			auto zz = atof(args[4].c_str())*units;
			
			insert(symbol, {xx, yy, zz});
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
		
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in the ions command" << std::endl;
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
