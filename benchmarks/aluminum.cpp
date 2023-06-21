/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

#include <sstream>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);

	int pardomains = 1;
	bool groundstate_only = false;
	vector3<int, contravariant> reps{2, 2, 2};
	int niter = 10;
	
	auto functional = input::interaction::pbe();
	
	{
		int opt;
		while ((opt = getopt(argc, argv, "p:?gs:yi:")) != EOF){
			switch(opt){
			case 'p':
				pardomains = atoi(optarg);
				break;
			case 'g':
				groundstate_only = true;
				break;
			case 's':
				sscanf(optarg, "%d%*c%d%*c%d", &reps[0], &reps[1], &reps[2]);
				break;
			case 'y':
				functional = input::interaction::pbe0();
				break;
 		  case 'i':
				niter = atoi(optarg);
				break;
			case '?':
				std::cerr << "usage is " << std::endl;
				std::cerr << "-p N to set the number of processors in the domain partition (1 by default)." << std::endl;
				std::cerr << "-g only calculate the ground state." << std::endl;
				std::cerr << "-s NxMxK the supercell size in each direction (2x2x2 by default)." << std::endl;
				std::cerr << "-y use the PBE0 hybrid functional (the default is PBE)." << std::endl;
				std::cerr << "-i N to set the number of SCF iterations to run (10 by default)." << std::endl;	
				exit(1);
			default:
				abort();
			}
		}
	}
	
	auto alat = 7.6524459_bohr;

	using frac_coord = vector3<double, contravariant>;
	
	std::vector<frac_coord> cell;

	cell.emplace_back(frac_coord{0.0, 0.0, 0.0});
	cell.emplace_back(frac_coord{0.0, 0.5, 0.5});
	cell.emplace_back(frac_coord{0.5, 0.0, 0.5});
	cell.emplace_back(frac_coord{0.5, 0.5, 0.0});
	
	systems::box box = systems::box::orthorhombic(reps[0]*alat, reps[1]*alat, reps[2]*alat);
	
	systems::ions ions(box);
	
	for(int ix = 0; ix < reps[0]; ix++){
		for(int iy = 0; iy < reps[1]; iy++){
			for(int iz = 0; iz < reps[2]; iz++){
				auto base = frac_coord{double(ix), double(iy), double(iz)};
				for(unsigned iatom = 0; iatom < cell.size(); iatom++){
					ions.insert_fractional("Al", (base + cell[iatom])/reps);
				}
			}
		}
	}
	
	assert(int(ions.geo().num_atoms()) == int(cell.size()*product(reps)));
				 
	systems::electrons electrons(env.par().states().domains(pardomains), ions, box, input::config::spacing(alat/20) | input::config::extra_states(2*product(reps)) | input::config::temperature(1000.0_K));
	
	auto restart_dir = "aluminum_" + std::to_string(reps[0]) + "_" + std::to_string(reps[1]) + "_" + std::to_string(reps[2]);

	auto not_found_gs = groundstate_only or not electrons.try_load(restart_dir);

	if(not_found_gs){
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, functional, inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(niter) | inq::input::scf::mixing(0.1));
		if(not groundstate_only) electrons.save(restart_dir);
	}

	if(not groundstate_only){
		real_time::propagate(ions, electrons, [](auto){}, functional, input::rt::num_steps(100) | input::rt::dt(0.05_atomictime));
	}
	
}

