/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);
	
	utils::match match(3.0e-5);

	auto box = systems::box::cubic(15.0_b).finite().cutoff_energy(30.0_Ha);

	systems::ions ions(box);

	auto distance = 121.0_pm;
	
	ions.insert("O", {-distance/2, 0.0_b, 0.0_b});
	ions.insert("O", {distance/2, 0.0_b, 0.0_b});	

	systems::electrons electrons(env.par(), ions, box, input::config::spin_polarized() | input::config::temperature(1000.0_K) | input::config::extra_states(3));
	ground_state::initial_guess(ions, electrons);
		
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), input::scf::mixing(0.1) | input::scf::energy_tolerance(1e-8_Ha));
	
	match.check("total energy",        result.energy.total()    ,  -32.983194928900);
	match.check("kinetic energy",      result.energy.kinetic()   ,  21.076512293160);
	match.check("eigenvalues",         result.energy.eigenvalues(), -7.232750608632);
	match.check("hartree energy",      result.energy.hartree(),     42.250517388010);
	match.check("external energy",     result.energy.external(),   -99.280824027265);
	match.check("non-local energy",    result.energy.nonlocal(),    -4.566033144402);
	match.check("XC energy",           result.energy.xc(),          -8.207482804083);
	match.check("XC density integral", result.energy.nvxc(),        -8.963440506145);
	match.check("HF exchange energy",  result.energy.hf_exchange(),  0.000000000000);
	match.check("ion-ion energy",      result.energy.ion(),         15.744115365679);

	return match.fail();
	
}

