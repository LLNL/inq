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
		
	electrons.load("al4h1_restart");

	std::vector<double> energy;
	
	real_time::propagate<>(ions, electrons, [&](auto data){energy.push_back(data.energy().total());}, input::interaction{}.lda(), options::real_time{}.num_steps(30).dt(0.055_atomictime));

	energy_match.check("energy step   0", energy[0],   -9.798687545996);
	energy_match.check("energy step  10", energy[10],  -9.798687876007);
	energy_match.check("energy step  20", energy[20],  -9.798688613754);
	energy_match.check("energy step  30", energy[30],  -9.798688812739);
	
	return energy_match.fail();
	
}

