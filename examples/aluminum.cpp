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

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <input/parse_xyz.hpp>
#include <config/path.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);

	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
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
	
	systems::ions ions(input::cell::orthorhombic(repx*alat*1.0_b, repy*alat*1.0_b, repz*alat*1.0_b), supercell);
	
	input::config conf;
	
	conf.extra_states = 8;
	conf.temperature = 300.0_K;
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0_Ha), conf);
	
	auto restart_dir = "aluminum_" + std::to_string(repx) + "_" + std::to_string(repy) + "_" + std::to_string(repz);
	
	ground_state::initialize(ions, electrons);

	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(200));
	
	energy_match.check("total energy",        result.energy.total()/supercell.size(),        -2.255100931633);
	energy_match.check("kinetic energy",      result.energy.kinetic()/supercell.size(),       0.917466797668);
	energy_match.check("eigenvalues",         result.energy.eigenvalues/supercell.size(),     0.246024442379);
	energy_match.check("Hartree energy",      result.energy.hartree/supercell.size(),         0.003056986010);
	energy_match.check("external energy",     result.energy.external/supercell.size(),        0.030210105865);
	energy_match.check("non-local energy",    result.energy.nonlocal/supercell.size(),        0.374607223743);
	energy_match.check("XC energy",           result.energy.xc/supercell.size(),             -1.081496717773);
	energy_match.check("XC density integral", result.energy.nvxc/supercell.size(),           -1.082373656917);
	energy_match.check("ion-ion energy",      result.energy.ion/supercell.size(),            -2.498945327145);

	inq::operations::io::save(restart_dir, electrons.phi_);

	return energy_match.fail();
	
}

