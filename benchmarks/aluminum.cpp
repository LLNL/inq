/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019-2022 Xavier Andrade, Alfredo A. Correa

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

#include <sstream>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);

	int pardomains = 1;
	bool groundstate_only = false;
	math::vector3<int> reps{2, 2, 2};

	auto functional = input::interaction::pbe();
	
	{
		int opt;
		while ((opt = getopt(argc, argv, "p:?gs:y")) != EOF){
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
			case '?':
				std::cerr << "usage is " << std::endl;
				std::cerr << "-p N to set the number of processors in the domain partition (1 by default)." << std::endl;
				std::cerr << "-g only calculate the ground state." << std::endl;
				std::cerr << "-s NxMxK the supercell size in each direction (2x2x2 by default)." << std::endl;
				std::cerr << "-y use the PBE0 hybrid functional (the default is PBE)." << std::endl;				
				exit(1);
			default:
				abort();
			}
		}
	}
	
	auto alat = 7.6524459_bohr;

	using frac_coord = math::vector3<decltype(0.0_crys)>;
	
	std::vector<frac_coord> cell;

	cell.emplace_back(frac_coord{0.0_crys, 0.0_crys, 0.0_crys});
	cell.emplace_back(frac_coord{0.0_crys, 0.5_crys, 0.5_crys});
	cell.emplace_back(frac_coord{0.5_crys, 0.0_crys, 0.5_crys});
	cell.emplace_back(frac_coord{0.5_crys, 0.5_crys, 0.0_crys});
	
	systems::box box = systems::box::orthorhombic(reps[0]*alat, reps[1]*alat, reps[2]*alat).spacing(alat/20);
	
	systems::ions ions(box);
	
	for(int ix = 0; ix < reps[0]; ix++){
		for(int iy = 0; iy < reps[1]; iy++){
			for(int iz = 0; iz < reps[2]; iz++){
				frac_coord base{ix*1.0_crys, iy*1.0_crys, iz*1.0_crys};
				for(unsigned iatom = 0; iatom < cell.size(); iatom++){
					ions.insert("Al", (base + cell[iatom])/reps);
				}
			}
		}
	}
	
	assert(int(ions.geo().num_atoms()) == int(cell.size()*product(reps)));
				 
	input::config conf;
				 
	conf.extra_states = 2*product(reps);
	conf.temperature = 300.0_K;	
	
	systems::electrons electrons(env.par().states().domains(pardomains), ions, box, conf);
	
	auto restart_dir = "aluminum_" + std::to_string(reps[0]) + "_" + std::to_string(reps[1]) + "_" + std::to_string(reps[2]);

	auto not_found_gs = groundstate_only or not electrons.load(restart_dir);

	if(not_found_gs){
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, functional, inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(10));
		electrons.save(restart_dir);
	}

	if(not groundstate_only){
		real_time::propagate(ions, electrons, [](auto){}, functional, input::rt::num_steps(100) | input::rt::dt(0.0565_atomictime));
	}
	
}

