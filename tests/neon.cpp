/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	auto comm = input::environment::global().comm();	

	utils::match energy_match(2.0e-5);

	{
		systems::ions ions(systems::cell::cubic(15.0_b).finite());
		ions.save(comm, ".default_ions");
	}

	{
		auto ions = systems::ions::load(".default_ions");		
		ions.insert(input::species("Ne"), {0.0_b, 0.0_b, 0.0_b});
		ions.save(comm, ".default_ions");
	}

	{
		auto el_opts = options::electrons{}.extra_states(3);
		el_opts.save(comm, ".default_electrons_options");
	}
	
	{
		auto el_opts = options::electrons::load(".default_electrons_options").cutoff(30.0_Ha);
		el_opts.save(comm, ".default_electrons_options");
	}

	{
		auto theo = options::theory{}.non_interacting();
		theo.save(comm, ".default_theory");
	}

	//REAL SPACE PSEUDO
	{
		auto ions = systems::ions::load(".default_ions");
		systems::electrons electrons(ions, options::electrons::load(".default_electrons_options"));
		
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, options::theory::load(".default_theory"));
		
		energy_match.check("total energy",     result.energy.total()      , -61.861045337100);
		energy_match.check("kinetic energy",   result.energy.kinetic()    ,  35.765610219604);
		energy_match.check("eigenvalues",      result.energy.eigenvalues(), -61.861045337100);
		energy_match.check("external energy",  result.energy.external()   , -79.509954154661);
		energy_match.check("non-local energy", result.energy.nonlocal()   , -18.116701402044);
		energy_match.check("ion-ion energy",   result.energy.ion()        ,   0.000000000000);

		electrons.save(".default_orbitals");
	}

	//FOURIER SPACE PSEUDO
	{
		auto el_opts = options::electrons::load(".default_electrons_options").fourier_pseudo();
		el_opts.save(comm, ".default_electrons_options");
	}
	
	{
		auto ions = systems::ions::load(".default_ions");
		systems::electrons electrons(ions, options::electrons::load(".default_electrons_options"));
		electrons.load(".default_orbitals");

		auto result = ground_state::calculate(ions, electrons, options::theory::load(".default_theory"));
		
		energy_match.check("total energy",     result.energy.total()      , -61.861056649453);
		energy_match.check("kinetic energy",   result.energy.kinetic()    ,  35.765555684056);
		energy_match.check("eigenvalues",      result.energy.eigenvalues(), -61.861056649453);
		energy_match.check("external energy",  result.energy.external()   , -79.509918897873);
		energy_match.check("non-local energy", result.energy.nonlocal()   , -18.116693435635);
		energy_match.check("ion-ion energy",   result.energy.ion()        ,   0.000000000000);
		
	}

	return energy_match.fail();
	
}
