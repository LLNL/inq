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
	
	utils::match energy_match(4.0e-6);

	double alat = 7.6524459;
	
	std::vector<input::atom> cell;

	cell.push_back( "Al" | alat*math::vector3<double>(0.0, 0.0, 0.0));
	cell.push_back( "Al" | alat*math::vector3<double>(0.0, 0.5, 0.5));
	cell.push_back( "Al" | alat*math::vector3<double>(0.5, 0.0, 0.5));
	cell.push_back( "Al" | alat*math::vector3<double>(0.5, 0.5, 0.0));	

	int repx = 2;
	int repy = 2;
	int repz = 2;

	std::vector<input::atom> supercell;

	for(int ix = 0; ix < repx; ix++){
		for(int iy = 0; iy < repy; iy++){
			for(int iz = 0; iz < repz; iz++){
				math::vector3<double> base{ix*alat, iy*alat, iz*alat};
				for(unsigned iatom = 0; iatom < cell.size(); iatom++){
					supercell.push_back(cell[iatom].species() | (base + cell[iatom].position()));
				}
			}
		}
	}

	assert(supercell.size() == cell.size()*repx*repy*repz);

	systems::box box = systems::box::orthorhombic(repx*alat*1.0_b, repy*alat*1.0_b, repz*alat*1.0_b).cutoff_energy(30.0_Ha);
	
	systems::ions ions(box, supercell);
	
	input::config conf;

	conf.extra_states = 2*repx*repy*repz;
	conf.temperature = 300.0_K;	
	
	systems::electrons electrons(env.dist().domains(1), ions, box, conf);
	
	auto restart_dir = "aluminum_" + std::to_string(repx) + "_" + std::to_string(repy) + "_" + std::to_string(repz);

	auto found_gs = electrons.load(restart_dir);

	if(not found_gs){

		// the parallelization distribution is different for the ground state
		systems::electrons gs_electrons(env.dist(), ions, box, conf);
		
		ground_state::initial_guess(ions, gs_electrons);
		
		auto result = ground_state::calculate(ions, gs_electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(200));
		
		gs_electrons.save(restart_dir);
		electrons.load(restart_dir);

	}

	auto propagation = real_time::propagate(ions, electrons, input::interaction::pbe(), input::rt::num_steps(100) | input::rt::dt(0.0625_atomictime));
	
	return energy_match.fail();
	
}

