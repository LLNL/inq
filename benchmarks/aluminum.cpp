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
	
	{
		int opt;
		while ((opt = getopt(argc, argv, "p:?")) != EOF){
			switch(opt){
			case 'p':
				pardomains = atoi(optarg);
				break;
			case '?':
				std::cerr << "usage is \n -p: <value> to set the number of processors in the domain partition." << std::endl;
				exit(0);
			default:
				abort();
			}
		}
	}
	
	utils::match energy_match(4.0e-6);

	auto alat = 7.6524459_bohr;

	using frac_coord = math::vector3<decltype(0.0_crys)>;
	
	std::vector<frac_coord> cell;

	cell.emplace_back(frac_coord{0.0_crys, 0.0_crys, 0.0_crys});
	cell.emplace_back(frac_coord{0.0_crys, 0.5_crys, 0.5_crys});
	cell.emplace_back(frac_coord{0.5_crys, 0.0_crys, 0.5_crys});
	cell.emplace_back(frac_coord{0.5_crys, 0.5_crys, 0.0_crys});

	int repx = 4;
	int repy = 4;
	int repz = 4;

	auto reps = math::vector3{repx, repy, repz};
	
	systems::box box = systems::box::orthorhombic(repx*alat, repy*alat, repz*alat).spacing(alat/20);
	
	systems::ions ions(box);
	
	for(int ix = 0; ix < repx; ix++){
		for(int iy = 0; iy < repy; iy++){
			for(int iz = 0; iz < repz; iz++){
				frac_coord base{ix*1.0_crys, iy*1.0_crys, iz*1.0_crys};
				for(unsigned iatom = 0; iatom < cell.size(); iatom++){
					ions.insert("Al", (base + cell[iatom])/reps);
				}
			}
		}
	}

	assert(int(ions.geo().num_atoms()) == int(cell.size()*repx*repy*repz));

	input::config conf;

	conf.extra_states = 2*repx*repy*repz;
	conf.temperature = 300.0_K;	
	
	systems::electrons electrons(env.par().states().domains(pardomains), ions, box, conf);
	
	auto restart_dir = "aluminum_" + std::to_string(repx) + "_" + std::to_string(repy) + "_" + std::to_string(repz);

	auto found_gs = electrons.load(restart_dir);

	if(not found_gs){
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(200));
	}

	auto propagation = real_time::propagate(ions, electrons, input::interaction::pbe(), input::rt::num_steps(100) | input::rt::dt(0.0565_atomictime));
	
	return energy_match.fail();
	
}

