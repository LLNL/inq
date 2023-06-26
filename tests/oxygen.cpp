/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);
	
	utils::match match(3.0e-5);

	systems::ions ions(ions::unit_cell::cubic(10.0_b).finite());

	auto distance = 121.0_pm;
	
	ions.insert("O", {-distance/2, 0.0_b, 0.0_b});
	ions.insert("O", {distance/2, 0.0_b, 0.0_b});	

	systems::electrons electrons(env.par(), ions, options::electrons{}.spacing(0.43_b).spin_polarized().temperature(1000.0_K).extra_states(2));
	ground_state::initial_guess(ions, electrons);
		
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), options::ground_state{}.mixing(0.2).energy_tolerance(1e-8_Ha));
	
	match.check("total energy",        result.energy.total()    ,   -32.885878270495);
	match.check("kinetic energy",      result.energy.kinetic()   ,   20.663840649123);
	match.check("eigenvalues",         result.energy.eigenvalues(),  -7.271345986319);
	match.check("hartree energy",      result.energy.hartree(),      42.107413515539);
	match.check("external energy",     result.energy.external(),    -98.967748472806);
	match.check("non-local energy",    result.energy.nonlocal(),     -4.258411086657);
	match.check("XC energy",           result.energy.xc(),           -8.175088241373);
	match.check("XC density integral", result.energy.nvxc(),         -8.923854107057);
	match.check("HF exchange energy",  result.energy.hf_exchange(),   0.000000000000);
	match.check("ion-ion energy",      result.energy.ion(),          15.744115365679);

	return match.fail();
	
}

