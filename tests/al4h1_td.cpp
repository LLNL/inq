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

	systems::box box = systems::box::cubic(alat).cutoff_energy(30.0_Ha);

	systems::ions ions(box);

	ions.insert("Al", {0.0_crys, 0.0_crys, 0.0_crys});
	ions.insert("Al", {0.0_crys, 0.5_crys, 0.5_crys});
	ions.insert("Al", {0.5_crys, 0.0_crys, 0.5_crys});
	ions.insert("Al", {0.5_crys, 0.5_crys, 0.0_crys});	
	ions.insert("H",  {0.1_crys, 0.2_crys, 0.3_crys});
	
	systems::electrons electrons(env.par(), ions, box,  input::config::extra_states(1) | input::config::temperature(300.0_K), input::kpoints::grid({2, 2, 2}, true));
	
	electrons.load("al4h1_restart");

	std::vector<double> energy;
	
	real_time::propagate<>(ions, electrons, [&](auto data){energy.push_back(data.energy());}, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));

	energy_match.check("energy step   0", energy[0],   -9.798687545996);
	energy_match.check("energy step  10", energy[10],  -9.798687876007);
	energy_match.check("energy step  20", energy[20],  -9.798688613754);
	energy_match.check("energy step  30", energy[30],  -9.798688812739);
	
	return energy_match.fail();
	
}

