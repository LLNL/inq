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

	std::vector<input::atom> geo;

	auto a = 10.18_b;

	auto box = systems::box::cubic(a);
	
	systems::ions ions(box);
	
	ions.insert_fractional("Si", {0.0,  0.0,  0.0 });
	ions.insert_fractional("Si", {0.25, 0.25, 0.25});
	ions.insert_fractional("Si", {0.5,  0.5,  0.0 });
	ions.insert_fractional("Si", {0.75, 0.75, 0.25});
	ions.insert_fractional("Si", {0.5,  0.0,  0.5 });
	ions.insert_fractional("Si", {0.75, 0.25, 0.75});
	ions.insert_fractional("Si", {0.0,  0.5,  0.5 });
	ions.insert_fractional("Si", {0.25, 0.75, 0.75});

	systems::electrons electrons(env.par(), ions, input::config::cutoff(35.0_Ha), input::kpoints::grid({1, 1, 1}, true));
	
	ground_state::initial_guess(ions, electrons);
	auto result = ground_state::calculate(ions, electrons, input::interaction::lda(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-9_Ha));
	electrons.save("silicon_restart");

	{ //NO KICK
		std::vector<double> jz;
		std::vector<double> Az;
		auto output = [&](auto data){
			jz.push_back(data.current()[2]); 
			Az.push_back(data.uniform_vector_potential()[2]);
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(40) | input::rt::dt(0.03_atomictime), ions::propagator::fixed{});
		
		data_match.check("current in z step   0", jz[0],   -0.000002678775);
		data_match.check("current in z step  10", jz[10],  -0.000002669488);
		data_match.check("current in z step  20", jz[20],  -0.000002663568);
		data_match.check("current in z step  30", jz[30],  -0.000002656289);
		data_match.check("current in z step  40", jz[40],  -0.000002637812);
		
		data_match.check("vector potential in z step   0", Az[0],   0.0);
		data_match.check("vector potential in z step  10", Az[10],  0.0);
		data_match.check("vector potential in z step  20", Az[20],  0.0);
		data_match.check("vector potential in z step  30", Az[30],  0.0);
		data_match.check("vector potential in z step  40", Az[40],  0.0);
	}
	
	{ //NO LRC CORRECTION
		electrons.load("silicon_restart");
				
		std::vector<double> jz;
		std::vector<double> Az;
		auto output = [&](auto data){
			jz.push_back(data.current()[2]); 
			Az.push_back(data.uniform_vector_potential()[2]);
		};
		
		auto kick = perturbations::kick{box.cell(), {0.0, 0.0, -0.005}, perturbations::gauge::velocity};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(40) | input::rt::dt(0.03_atomictime), ions::propagator::fixed{}, kick);
		
		data_match.check("current in z step   0", jz[0],   -0.157729547895);
		data_match.check("current in z step  10", jz[10],  -0.151911067747);
		data_match.check("current in z step  20", jz[20],  -0.144088608186);
		data_match.check("current in z step  30", jz[30],  -0.140191267000);
		data_match.check("current in z step  40", jz[40],  -0.135683835952);
		
		data_match.check("vector potential in z step   0", Az[0],   0.050900000000);
		data_match.check("vector potential in z step  10", Az[10],  0.050900000000);
		data_match.check("vector potential in z step  20", Az[20],  0.050900000000);
		data_match.check("vector potential in z step  30", Az[30],  0.050900000000);
		data_match.check("vector potential in z step  40", Az[40],  0.050900000000);
	}

	{ //LRC CORRECTION
		electrons.load("silicon_restart");
		
		std::vector<double> jz;
		std::vector<double> Az;
		auto output = [&](auto data){
			jz.push_back(data.current()[2]); 
			Az.push_back(data.uniform_vector_potential()[2]);
		};
		
		auto kick = perturbations::kick{box.cell(), {0.0, 0.0, -0.005}, perturbations::gauge::velocity};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda() | input::interaction::lrc(0.2), input::rt::num_steps(40) | input::rt::dt(0.03_atomictime), ions::propagator::fixed{}, kick);
		
		data_match.check("current in z step   0", jz[0],   -0.157729547895);
		data_match.check("current in z step  10", jz[10],  -0.151948847224);
		data_match.check("current in z step  20", jz[20],  -0.144242611185);
		data_match.check("current in z step  30", jz[30],  -0.140532611854);
		data_match.check("current in z step  40", jz[40],  -0.136279338368);
		
		data_match.check("vector potential in z step   0", Az[0],   0.050900000000);
		data_match.check("vector potential in z step  10", Az[10],  0.050912284789);
		data_match.check("vector potential in z step  20", Az[20],  0.050950852243);
		data_match.check("vector potential in z step  30", Az[30],  0.051014538212);
		data_match.check("vector potential in z step  40", Az[40],  0.051102615912);
	}
	
	fftw_cleanup();

	return data_match.fail();
}

