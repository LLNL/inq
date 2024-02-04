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
	
	utils::match energy_match(2.0e-5);

	auto cell = systems::cell::cubic(15.0_b).finite();

	{
		systems::ions ions(cell);
		
		ions.insert(input::species("C").pseudo(inq::config::path::unit_tests_data() + "C_ONCV_PBE-1.2.xml"), {0.0_b, 0.0_b, 0.0_b});

		systems::electrons electrons(ions, options::electrons{}.cutoff(25.0_Ha).extra_states(4).temperature(300.0_K));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total()    ,   -5.364150140372);
		energy_match.check("kinetic energy",      result.energy.kinetic()   ,   3.176958479664);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(), -1.404891046908);
		energy_match.check("Hartree energy",      result.energy.hartree(),      4.371229020227);
		energy_match.check("external energy",     result.energy.external(),   -11.993822927726);
		energy_match.check("non-local energy",    result.energy.nonlocal(),     0.494166478064);
		energy_match.check("XC energy",           result.energy.xc(),          -1.412681190601);
		energy_match.check("XC density integral", result.energy.nvxc(),        -1.824651117364);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),          0.000000000000);

	}

	{
		systems::ions ions(cell);
		
		ions.insert(input::species("C"), {0.0_b, 0.0_b, 0.0_b});

		systems::electrons electrons(ions, options::electrons{}.pseudopotentials(pseudo::set::ccecp()).cutoff(25.0_Ha).extra_states(4).temperature(300.0_K));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total()    ,   -5.390943623548);
		energy_match.check("kinetic energy",      result.energy.kinetic()  ,    3.462168739143);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(), -1.409268764339);
		energy_match.check("Hartree energy",      result.energy.hartree(),      4.396047239823);
		energy_match.check("external energy",     result.energy.external(),   -12.434345771131);
		energy_match.check("non-local energy",    result.energy.nonlocal(),     0.607257611969);
		energy_match.check("XC energy",           result.energy.xc(),          -1.422071443352);
		energy_match.check("XC density integral", result.energy.nvxc(),        -1.836443823966);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),          0.000000000000);
	}
	return energy_match.fail();
	
}
