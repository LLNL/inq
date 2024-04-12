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
	
	constexpr auto name() const {
		return "perturbations";
	}

	constexpr auto one_line() const {
		return "Defines the perturbations in a real-time simulation";
	}
	constexpr auto help() const {
		
		return R""""(

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


- `perturbations laser <px> <py> <pz> <frequency|wavelength> <f> <units>`

  Adds a laser perturbation. The laser is approximated by a
  monochromatic electric field (the magnetic part is ignored). The
  direction and intensity of the field is given by the values _px_,
  _py_ and _pz_.

  You also have to specify the frequency of the laser using the
  keyword 'frequency' followed by its value and its units (you can use
  'eV' or 'THz' for example, check `inq help units` for a full list of
  the energy units available).

  Alternative you can use the keyword 'wavelength' followed by the
  wavelength of the laser and its units (typically 'nm').

  For finite systems the laser is applied in the length gauge as a
  potential. For periodic systems, the laser is applied in the
  velocity gauge as a uniform vector potential.

  Examples: `inq perturbations laser  0    0    0.5 frequency  300 THz`
            `inq perturbations laser -0.01 0.01 0   wavelength 880 nm`


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

	void laser(vector3<double, cartesian> const & polarization, quantity<magnitude::energy> freq) const {
		auto perturbations = perturbations::blend::load(".inq/default_perturbations");
		auto cell = systems::ions::load(".inq/default_ions").cell();
		if(cell.periodicity() == 0){
			perturbations.add(perturbations::laser(polarization, freq, perturbations::gauge::length));
		} else {
			perturbations.add(perturbations::laser(polarization, freq, perturbations::gauge::velocity));
		}
		perturbations.save(input::environment::global().comm(), ".inq/default_perturbations");
	}
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		using namespace magnitude;
		using utils::str_to;
		
		if(args.size() == 0) {
			operator()();
			actions::normal_exit();
		}

		if(args[0] == "clear"){
			if(args.size() != 1) actions::error(input::environment::global().comm(), "The 'perturbations clear' command doesn't take arguments");
			clear();
			if(not quiet) operator()();
			actions::normal_exit();
		}

		if(args[0] == "kick"){
			if(args.size() != 4) actions::error(input::environment::global().comm(), "Invalid arguments for 'perturbations kick'");
			kick({str_to<double>(args[1]), str_to<double>(args[2]), str_to<double>(args[3])});
			if(not quiet) operator()();
			actions::normal_exit();
		}

		if(args[0] == "laser"){
			if(args.size() != 7) actions::error(input::environment::global().comm(), "Invalid syntax for 'perturbations laser'.\nUse 'inq perturbations laser <px> <py> <pz> frequency <f> <units>'");

			quantity<energy> freq;
			
			if(args[4] == "frequency") {
				freq = energy::parse(str_to<double>(args[5]), args[6]);
			} else if(args[4] == "wavelength") {
				freq = energy::from_wavelength(length::parse(str_to<double>(args[5]), args[6]));
			} else {
				actions::error(input::environment::global().comm(), "Invalid argument " + args[4] + " for 'perturbations laser'.\nUse 'inq perturbations laser <px> <py> <pz> <frequency|wavelength> <f> <units>'");
			}

			laser({str_to<double>(args[1]), str_to<double>(args[2]), str_to<double>(args[3])}, freq);
			if(not quiet) operator()();
			actions::normal_exit();
		}

		actions::error(input::environment::global().comm(), "Invalid syntax in the perturbations command");
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
