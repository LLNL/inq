/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <observables/kinetic_energy_density.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);

	utils::match energy_match(3.0e-5);

	auto dcc = 1.42_A;
  auto aa = sqrt(3)*dcc;
  auto lz = 10.0_b;

	systems::ions ions(ions::unit_cell::lattice(aa*vector3{1.0, 0.0, 0.0}, aa*vector3{-1.0/2.0, sqrt(3.0)/2.0, 0.0}, {0.0_b, 0.0_b, lz}).periodicity(2));
	
	ions.insert("C", {0.0_b, 0.0_b, 0.0_b});
	ions.insert("C", {0.0_b, dcc,   0.0_b});

	{
		systems::electrons electrons(env.par(), ions, options::electrons{}.spacing(aa/15.0).extra_states(2), input::kpoints::grid({1, 1, 1}, false));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -11.794106663282);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         9.555702987386);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -3.535844485436);
		energy_match.check("Hartree energy",      result.energy.hartree(),       -11.139895259643);
		energy_match.check("external energy",     result.energy.external(),       15.119948124584);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.139268347279);
		energy_match.check("XC energy",           result.energy.xc(),             -4.376877703296);
		energy_match.check("XC density integral", result.energy.nvxc(),           -4.792436730841);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -19.813716465033);
		
	}
		
	fftw_cleanup();
	
	return energy_match.fail();
	
}


