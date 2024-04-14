/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__IONS
#define INQ__INTERFACE__IONS

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <interface/actions.hpp>
#include <interface/cell.hpp>
#include <systems/ions.hpp>

#ifdef INQ_PYTHON_INTERFACE
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

namespace inq {
namespace interface {

struct {
	
	constexpr auto name() const {
		return "ions";
	}

	constexpr auto one_line() const {
		return "Defines the ions in the simulation";
	}
	constexpr auto help() const {
		
		return R""""(

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

	static void show() {
		auto ions = systems::ions::load(".inq/default_ions");		
		if(input::environment::global().comm().root()) std::cout << ions;
	}

	void operator()() const {
		show();
	}

	static void insert(input::species const & sp, vector3<quantity<magnitude::length>> const & pos) {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.insert(sp, pos);
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	static void insert_fractional(input::species const & sp, vector3<double, contravariant> const & pos) {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.insert_fractional(sp, pos);
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	static void clear() {
		auto ions = systems::ions::load(".inq/default_ions");
		ions.clear();
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	static void file(std::string const & filename) {
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
	
	static void file(std::string const & filename, quantity<magnitude::length> const & radius) {
		auto ions = systems::ions::parse(filename, radius);
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

		using utils::str_to;
		
		if(args.size() == 0) {
			operator()();
			actions::normal_exit();
		}

		if(args[0] == "clear"){

			if(args.size() != 1) actions::error(input::environment::global().comm(), "The 'ions clear' command doesn't take arguments.");
			clear();
			if(not quiet) operator()();
			actions::normal_exit();
		}

 		if(args.size() >= 2 and args[0] == "insert" and args[1] == "fractional"){

			if(args.size() != 6) actions::error(input::environment::global().comm(), "Wrong arguments for ions insert.\nUse: inq ions insert fractional <symbol> <x> <y> <z>");
			
			auto symbol = args[2];
			auto xx = str_to<double>(args[3]);
			auto yy = str_to<double>(args[4]);
			auto zz = str_to<double>(args[5]);
			
			insert_fractional(symbol, {xx, yy, zz});
			if(not quiet) operator()();
			actions::normal_exit();
		}

		if(args[0] == "insert"){

			if(args.size() != 6) actions::error(input::environment::global().comm(), "Wrong arguments for ions insert.\nUse: inq ions insert <symbol> <x> <y> <z> <units>");

			auto symbol = args[1];
			auto units = args[5];
			auto xx = magnitude::length::parse(str_to<double>(args[2]), units);
			auto yy = magnitude::length::parse(str_to<double>(args[3]), units);
			auto zz = magnitude::length::parse(str_to<double>(args[4]), units);
			
			insert(symbol, {xx, yy, zz});
			if(not quiet) operator()();
			actions::normal_exit();
		}

		if(args.size() == 2 and args[0] == "file"){
			file(args[1]);
			if(not quiet) {
				interface::cell();
				operator()();
			}
			actions::normal_exit();
		}
		
		if(args.size() == 5 and args[0] == "file" and args[2] == "radius"){
			auto radius = magnitude::length::parse(str_to<double>(args[3]), args[4]);
			
			file(args[1], radius);
			if(not quiet) {
				interface::cell();
				operator()();
			}
			actions::normal_exit();
		}
		
		if(input::environment::global().comm().root()) actions::error(input::environment::global().comm(), "Invalid syntax in the ions command");
	}

#ifdef INQ_PYTHON_INTERFACE
	template <class PythonModule>
	void python_interface(PythonModule & module) const {
		namespace py = pybind11;
		using namespace pybind11::literals;

		auto sub = module.def_submodule(name(), help());

		sub.def("show", &show);

		sub.def("clear", &clear);

		sub.def("insert", [](std::string const & symbol, std::vector<double> const & coords, std::string const & units) {
			auto c0 = magnitude::length::parse(coords[0], units);
			auto c1 = magnitude::length::parse(coords[1], units);
			auto c2 = magnitude::length::parse(coords[2], units);

			insert(symbol, {c0, c1, c2});
		}, "species"_a,  "coordinates"_a,  "units"_a);

		sub.def("insert_fractional", [](std::string const & symbol, std::vector<double> const & coords) {
			insert_fractional(symbol, {coords[0], coords[1], coords[2]});
		}, "species"_a,  "coordinates"_a);

		sub.def("file", [](std::string const & filename) {
			file(filename);
		}, "filename"_a);
		
		sub.def("file", [](std::string const & filename, double radius, std::string const & units) {
			file(filename, magnitude::length::parse(radius, units));
		}, "filename"_a, "radius"_a, "units"_a);
		
	}
#endif
	
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
