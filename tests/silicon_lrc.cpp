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

	utils::match data_match(7.0e-5);

	auto a = 10.18_b;
	systems::ions ions(systems::cell::cubic(a));
	
	ions.insert_fractional("Si", {0.0,  0.0,  0.0 });
	ions.insert_fractional("Si", {0.25, 0.25, 0.25});
	ions.insert_fractional("Si", {0.5,  0.5,  0.0 });
	ions.insert_fractional("Si", {0.75, 0.75, 0.25});
	ions.insert_fractional("Si", {0.5,  0.0,  0.5 });
	ions.insert_fractional("Si", {0.75, 0.25, 0.75});
	ions.insert_fractional("Si", {0.0,  0.5,  0.5 });
	ions.insert_fractional("Si", {0.25, 0.75, 0.75});

	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(35.0_Ha), input::kpoints::grid({1, 1, 1}, true));
	
	ground_state::initial_guess(ions, electrons);
	auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), options::ground_state{}.energy_tolerance(1e-9_Ha));

	data_match.check("total energy",        result.energy.total(),          -33.418896726864);
  data_match.check("kinetic energy",      result.energy.kinetic(),         13.130726145119);
  data_match.check("eigenvalues",         result.energy.eigenvalues(),     -0.213554205314);
  data_match.check("Hartree energy",      result.energy.hartree(),          2.381711162124);
  data_match.check("external energy",     result.energy.external(),        -9.463654074017);
  data_match.check("non-local energy",    result.energy.nonlocal(),         4.498591799053);
  data_match.check("XC energy",           result.energy.xc(),             -12.482651264043);
  data_match.check("XC density integral", result.energy.nvxc(),           -13.142640399717);
  data_match.check("ion-ion energy",      result.energy.ion(),            -31.483620495100);
	
	electrons.save("silicon_lrc_restart");

	{ //NO KICK
		std::vector<double> jz;
		std::vector<double> Az;
		std::vector<double> energy;
		auto output = [&](auto data){
			jz.push_back(data.current()[2]); 
			Az.push_back(data.uniform_vector_potential()[2]);
			energy.push_back(data.energy().total());
		};
		
		real_time::propagate<>(ions, electrons, output, options::theory{}.lda(), options::real_time{}.num_steps(40).dt(0.03_atomictime), ions::propagator::fixed{});

		data_match.check("energy step   0", energy[0],   -33.418896726864);
		data_match.check("energy step  10", energy[10],  -33.418896726864);
		data_match.check("energy step  20", energy[20],  -33.418896726864);
		data_match.check("energy step  30", energy[30],  -33.418896726864);
		data_match.check("energy step  40", energy[40],  -33.418896726864);
		
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
		electrons.load("silicon_lrc_restart");
				
		std::vector<double> jz;
		std::vector<double> Az;
		std::vector<double> energy;
		auto output = [&](auto data){
			jz.push_back(data.current()[2]); 
			Az.push_back(data.uniform_vector_potential()[2]);
			energy.push_back(data.energy().total());
		};
		
		auto kick = perturbations::kick{ions.cell(), {0.0, 0.0, -0.005}, perturbations::gauge::velocity};
		
		real_time::propagate<>(ions, electrons, output, options::theory{}.lda(), options::real_time{}.num_steps(40) .dt(0.03_atomictime), ions::propagator::fixed{}, kick);

		data_match.check("energy step   0", energy[0],   -33.418518667909);
		data_match.check("energy step  10", energy[10],  -33.418518667745);
		data_match.check("energy step  20", energy[20],  -33.418518667170);
		data_match.check("energy step  30", energy[30],  -33.418518666925);
		data_match.check("energy step  40", energy[40],  -33.418518667053);
		
		data_match.check("current in z step   0", jz[0],   -0.157728627231);
		data_match.check("current in z step  10", jz[10],  -0.151910113914);
		data_match.check("current in z step  20", jz[20],  -0.144087622581);
		data_match.check("current in z step  30", jz[30],  -0.140190241687);
		data_match.check("current in z step  40", jz[40],  -0.135682773098);
		
		data_match.check("vector potential in z step   0", Az[0],   0.050900000000);
		data_match.check("vector potential in z step  10", Az[10],  0.050900000000);
		data_match.check("vector potential in z step  20", Az[20],  0.050900000000);
		data_match.check("vector potential in z step  30", Az[30],  0.050900000000);
		data_match.check("vector potential in z step  40", Az[40],  0.050900000000);
	}

	{ //LRC CORRECTION
		electrons.load("silicon_lrc_restart");
		
		std::vector<double> jz;
		std::vector<double> Az;
		std::vector<double> energy;		
		auto output = [&](auto data){
			jz.push_back(data.current()[2]); 
			Az.push_back(data.uniform_vector_potential()[2]);
			energy.push_back(data.energy().total());
		};
		
		auto kick = perturbations::kick{ions.cell(), {0.0, 0.0, -0.005}, perturbations::gauge::velocity};
		
		real_time::propagate<>(ions, electrons, output, options::theory{}.lda().lrc(0.2), options::real_time{}.num_steps(40).dt(0.03_atomictime), ions::propagator::fixed{}, kick);

		data_match.check("energy step   0", energy[0],   -33.418518667909);
		data_match.check("energy step  10", energy[10],  -33.418518487977);
		data_match.check("energy step  20", energy[20],  -33.418517946804);
		data_match.check("energy step  30", energy[30],  -33.418517094998);
		data_match.check("energy step  40", energy[40],  -33.418515958018);
			
		data_match.check("current in z step   0", jz[0],   -0.157728627231);
		data_match.check("current in z step  10", jz[10],  -0.151948001118);
		data_match.check("current in z step  20", jz[20],  -0.144242243331);
		data_match.check("current in z step  30", jz[30],  -0.140532994159);
		data_match.check("current in z step  40", jz[40],  -0.136280648873);
		
		data_match.check("vector potential in z step   0", Az[0],   0.050900000000);
		data_match.check("vector potential in z step  10", Az[10],  0.050912319754);
		data_match.check("vector potential in z step  20", Az[20],  0.050951055407);
		data_match.check("vector potential in z step  30", Az[30],  0.051015009633);
		data_match.check("vector potential in z step  40", Az[40],  0.051103423993);
	}
	
	fftw_cleanup();

	return data_match.fail();
}

