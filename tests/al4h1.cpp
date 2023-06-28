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
	
	input::environment env(argc, argv);
	
	utils::match energy_match(1.0e-5);

	auto alat = 7.6524459_bohr;
	systems::ions ions(systems::cell::cubic(alat));

	ions.insert_fractional("Al", {0.0, 0.0, 0.0});
	ions.insert_fractional("Al", {0.0, 0.5, 0.5});
	ions.insert_fractional("Al", {0.5, 0.0, 0.5});
	ions.insert_fractional("Al", {0.5, 0.5, 0.0});	
	ions.insert_fractional("H",  {0.1, 0.2, 0.3});
	
	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(1).temperature(300.0_K), input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), options::ground_state{}.energy_tolerance(1e-8_Ha).calculate_forces());
	
	energy_match.check("total energy",        result.energy.total(),        -9.802337517464);
	energy_match.check("kinetic energy",      result.energy.kinetic(),       4.200428677371);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),   0.602437060724);
	energy_match.check("Hartree energy",      result.energy.hartree(),       0.219185000916);
	energy_match.check("external energy",     result.energy.external(),     -0.562801499284);
	energy_match.check("non-local energy",    result.energy.nonlocal(),      1.427219744575);
	energy_match.check("XC energy",           result.energy.xc(),           -4.767997327810);
	energy_match.check("XC density integral", result.energy.nvxc(),         -4.900779863769);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange(),   0.000000000000);
	energy_match.check("ion-ion energy",      result.energy.ion(),         -10.318372113231);
	
	energy_match.check("force 1 x",           result.forces[0][0],          -0.022482609447);
	energy_match.check("force 1 y",           result.forces[0][1],          -0.041215135710);
	energy_match.check("force 1 z",           result.forces[0][2],          -0.052722871506);
	energy_match.check("force 2 x",           result.forces[1][0],          -0.022475634954);
	energy_match.check("force 2 y",           result.forces[1][1],           0.052695862055);
	energy_match.check("force 2 z",           result.forces[1][2],           0.041206689231);
	energy_match.check("force 3 x",           result.forces[2][0],           0.005729211291);
	energy_match.check("force 3 y",           result.forces[2][1],          -0.012778639911);
	energy_match.check("force 3 z",           result.forces[2][2],           0.012775163241);
	energy_match.check("force 4 x",           result.forces[3][0],           0.007075771429);
	energy_match.check("force 4 y",           result.forces[3][1],           0.012276024136);
	energy_match.check("force 4 z",           result.forces[3][2],          -0.012279911145);
	energy_match.check("force 5 x",           result.forces[4][0],           0.027655896120);
	energy_match.check("force 5 y",           result.forces[4][1],          -0.010188638343);
	energy_match.check("force 5 z",           result.forces[4][2],           0.010361037452);
	
	electrons.save("al4h1_restart");
	
	return energy_match.fail();
	
}

