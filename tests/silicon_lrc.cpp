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

	utils::match data_match(3.0e-5);

	auto a = 10.18_b;
	systems::ions ions(ions::unit_cell::cubic(a));
	
	ions.insert_fractional("Si", {0.0,  0.0,  0.0 });
	ions.insert_fractional("Si", {0.25, 0.25, 0.25});
	ions.insert_fractional("Si", {0.5,  0.5,  0.0 });
	ions.insert_fractional("Si", {0.75, 0.75, 0.25});
	ions.insert_fractional("Si", {0.5,  0.0,  0.5 });
	ions.insert_fractional("Si", {0.75, 0.25, 0.75});
	ions.insert_fractional("Si", {0.0,  0.5,  0.5 });
	ions.insert_fractional("Si", {0.25, 0.75, 0.75});

	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(25.0_Ha), input::kpoints::grid({1, 1, 1}, true));
	
	ground_state::initial_guess(ions, electrons);

	auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-9_Ha));
	
	std::vector<double> jz;
	std::vector<double> Az;
	auto output = [&](auto data){
		jz.push_back(data.current()[2]); 
		Az.push_back(data.vector_field()[2]);
	};

	auto kick = perturbations::kick{ions.cell(), {0.0, 0.0, -0.005}, perturbations::gauge::velocity};

	real_time::propagate<>(ions, electrons, output, input::interaction::lda()|input::interaction::gauge_field(0.2), input::rt::num_steps(30) | input::rt::dt(0.04_atomictime), ions::propagator::fixed{}, kick);

	data_match.check("current in z step   0", jz[0],   -0.157798692934);
	data_match.check("current in z step  10", jz[10],  -0.146691801751);
	data_match.check("current in z step  20", jz[20],  -0.138123785595);
	data_match.check("current in z step  30", jz[30],  -0.132646541020);
	
	data_match.check("vector potential in z step   0", Az[0],   0.050900000000);
	data_match.check("vector potential in z step  10", Az[10],  0.050921665475);
	data_match.check("vector potential in z step  20", Az[20],  0.050988650385);
	data_match.check("vector potential in z step  30", Az[30],  0.051098335228);
	
	fftw_cleanup();

	return data_match.fail();
}

