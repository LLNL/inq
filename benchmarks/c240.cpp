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

	int pardomains = 1;
	bool groundstate_only = false;

	{
		int opt;
		while ((opt = getopt(argc, argv, "p:?gs:")) != EOF){
			switch(opt){
			case 'p':
				pardomains = atoi(optarg);
				break;
			case 'g':
				groundstate_only = true;
				break;
			case '?':
				std::cerr << "usage is " << std::endl;
				std::cerr << "-p N to set the number of processors in the domain partition (1 by default)." << std::endl;
				std::cerr << "-g only calculate the ground state." << std::endl;				
				exit(1);
			default:
				abort();
			}
		}
	}
	
	auto box = systems::box::cubic(20.0_A).finite();
	
	systems::ions ions(box);

	ions.insert(input::parse_xyz(config::path::unit_tests_data() + "c240.xyz"));
	
	std::string restart_dir = "c240_restart";
	
	systems::electrons electrons(env.par().states().domains(pardomains), ions, box, input::config::spacing(20.0_A/90) | input::config::extra_states(32) | input::config::temperature(300.0_K));

	auto not_found_gs = groundstate_only or not electrons.try_load(restart_dir);
		
	if(not_found_gs){
		ground_state::initial_guess(ions, electrons);
		ground_state::calculate(ions, electrons, input::interaction::pbe(), input::scf::steepest_descent() | input::scf::scf_steps(10) | input::scf::mixing(0.3));
		electrons.save(restart_dir);
	}
	
	if(not groundstate_only){
		real_time::propagate(ions, electrons, [](auto){}, input::interaction::pbe(), input::rt::num_steps(100) | input::rt::dt(0.0565_atomictime));
	}
	
}

