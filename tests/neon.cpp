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
		
	utils::match energy_match(2.0e-5);

	systems::ions ions(ions::unit_cell::cubic(15.0_b).finite());
	ions.insert("Ne" | input::species::nofilter(), {0.0_b, 0.0_b, 0.0_b});

	systems::electrons electrons(env.par(), ions, input::config::extra_states(3) | input::config::cutoff(30.0_Ha));
	
	ground_state::initial_guess(ions, electrons);
	
	//REAL SPACE PSEUDO
	{
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting());
		
		energy_match.check("total energy",     result.energy.total()      , -61.765371991220);
		energy_match.check("kinetic energy",   result.energy.kinetic()    ,  35.606535224997);
		energy_match.check("eigenvalues",      result.energy.eigenvalues(), -61.765371991153);
		energy_match.check("external energy",  result.energy.external()   , -79.444978433236);
		energy_match.check("non-local energy", result.energy.nonlocal()   , -17.926928782914);
		energy_match.check("ion-ion energy",   result.energy.ion()        ,   0.000000000000);

	}

	//FOURIER SPACE PSEUDO
	{	
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting() | input::interaction::fourier_pseudo());
		
		energy_match.check("total energy",     result.energy.total()      , -61.765376105880);
		energy_match.check("kinetic energy",   result.energy.kinetic()    ,  35.606511739929);
		energy_match.check("eigenvalues",      result.energy.eigenvalues(), -61.765376105880);
		energy_match.check("external energy",  result.energy.external()   , -79.444966794608);
		energy_match.check("non-local energy", result.energy.nonlocal()   , -17.926921051201);
		energy_match.check("ion-ion energy",   result.energy.ion()        ,   0.000000000000);
		
	}

	return energy_match.fail();
	
}
