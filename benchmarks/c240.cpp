/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

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
	
	auto box = systems::box::cubic(20.0_A).finite().spacing(20.0_A/90);
	
	systems::ions ions(box);

	ions.insert(input::parse_xyz(config::path::unit_tests_data() + "c240.xyz"));
	
	std::string restart_dir = "c240_restart";
	
	systems::electrons electrons(env.par().states().domains(pardomains), ions, box, input::config::extra_states(32) | input::config::temperature(300.0_K));

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

