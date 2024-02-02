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

	auto & env = inq::input::environment::global();

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
	
	auto ions = systems::ions::parse(inq::config::path::unit_tests_data() + "c240.xyz", systems::cell::cubic(20.0_A).finite());
		
	std::string restart_dir = "c240_restart";
	
	systems::electrons electrons(env.par().states().domains(pardomains), ions, options::electrons{}.spacing(20.0_A/90).extra_states(32).temperature(300.0_K));

	auto not_found_gs = groundstate_only or not electrons.try_load(restart_dir);
		
	if(not_found_gs){
		ground_state::initial_guess(ions, electrons);
		try { ground_state::calculate(ions, electrons, options::theory{}.pbe(), options::ground_state{}.steepest_descent().max_steps(10).mixing(0.3)); }
		catch(...){ }
		electrons.save(restart_dir);
	}
	
	if(not groundstate_only){
		real_time::propagate(ions, electrons, [](auto){}, options::theory{}.pbe(), options::real_time{}.num_steps(100).dt(0.0565_atomictime));
	}
	
}

