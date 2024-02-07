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

This command defines the kpoints to be used in the simulation. This
command has to be given after the ions are defined.

Note that right now inq doesn't use symmetries to reduce the kpoint
grid.

These are the options available:

- `kpoints`

   Prints a list of the current kpoints in the simulation.

   Example: `inq kpoints`


- `kpoints gamma` (default)

   Tells the system to only use the 0 0 0 point. This is the default.

   Example: `inq kpoints gamma`


- `kpoints grid <nx> <ny> <nz>`

   Define a uniform grid to sample the Brillouin zone with sizes _nx_,
   _ny_ and _nz_. This generates a non-shifted grid that will include
   gamma.

   Examples: `inq kpoints grid 8 8 8`


- `kpoints shifted grid <nx> <ny> <nz>`

   Define a shifted uniform grid to sample the Brillouin zone with sizes _nx_,
   _ny_ and _nz_. This grid doesn't include gamma.

   Examples: `inq kpoints shifted grid 4 4 4`


- `kpoints insert <kx> <ky> <kz> <w>`

   Add a kpoint to the grid with coordinates _kx_, _ky_, _kz_, and
   weight _w_. The kpoint must be given in covariant coordinates in the
   range [-0.5, 0.5).

   Examples: `inq kpoints insert 0.25 0.25 0.25 1.0`


- `kpoints clear`

   Removes all the kpoints defined in the grid.

   Examples: `inq kpoints clear`


)"""";
	}

  void operator()() const {
		auto bz = ions::brillouin::load(".inq/default_brillouin");
		if(input::environment::global().comm().root()) std::cout << bz;
	}

  void gamma() const {
		auto bz = ions::brillouin(systems::ions::load(".inq/default_ions"), input::kpoints::gamma());
		bz.save(input::environment::global().comm(), ".inq/default_brillouin");
	}

  void grid(int nx, int ny, int nz) const {
		auto bz = ions::brillouin(systems::ions::load(".inq/default_ions"), input::kpoints::grid({nx, ny, nz}));
		bz.save(input::environment::global().comm(), ".inq/default_brillouin");
	}

  void shifted_grid(int nx, int ny, int nz) const {
		auto bz = ions::brillouin(systems::ions::load(".inq/default_ions"), input::kpoints::grid({nx, ny, nz}, true));
		bz.save(input::environment::global().comm(), ".inq/default_brillouin");
	}

  void insert(double const & kx, double const & ky, double const & kz, double const & weight) const {
    auto bz = ions::brillouin{};
    try { bz = ions::brillouin::load(".inq/default_brillouin"); }
    catch(...){  }
    bz.insert({kx, ky, kz}, weight);
		bz.save(input::environment::global().comm(), ".inq/default_brillouin");
  }

  void clear() const {
    auto bz = ions::brillouin{};
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

    if(args.size() == 4 and args[0] == "grid") {
      
      grid(atoi(args[1].c_str()), atoi(args[2].c_str()), atoi(args[3].c_str()));
      if(not quiet) operator()();
			exit(0);
		}

    if(args.size() == 4 and args[0] == "shifted-grid") {
      
      shifted_grid(atoi(args[1].c_str()), atoi(args[2].c_str()), atoi(args[3].c_str()));
      if(not quiet) operator()();
			exit(0);
		}
    
    if(args.size() == 5 and args[0] == "insert") {
      
      insert(atof(args[1].c_str()), atof(args[2].c_str()), atof(args[3].c_str()), atof(args[4].c_str()));
      if(not quiet) operator()();
			exit(0);
		}
    
    if(args.size() == 1 and args[0] == "clear") {
			clear();
      if(not quiet) operator()();
			exit(0);
		}
    
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in the 'kpoints' command" << std::endl;
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
