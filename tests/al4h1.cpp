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

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {boost::mpi3::fill, 1});
	
	utils::match energy_match(4.0e-6);

	double alat = 7.6524459;
	
	std::vector<input::atom> cell;

	cell.push_back( "Al" | alat*math::vector3<double>(0.0, 0.0, 0.0));
	cell.push_back( "Al" | alat*math::vector3<double>(0.0, 0.5, 0.5));
	cell.push_back( "Al" | alat*math::vector3<double>(0.5, 0.0, 0.5));
	cell.push_back( "Al" | alat*math::vector3<double>(0.5, 0.5, 0.0));	
	cell.push_back( "H"  | alat*math::vector3<double>(0.1, 0.2, 0.3));
	
	systems::box box = systems::box::orthorhombic(alat*1.0_b, alat*1.0_b, alat*1.0_b).cutoff_energy(30.0_Ha);
	systems::ions ions(box, cell);
	
	input::config conf;

	conf.extra_states = 1;
	conf.temperature = 300.0_K;	
	
	systems::electrons electrons(cart_comm, ions, box, conf, input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe());

	energy_match.check("ion-ion energy",      result.energy.ion,             -10.318372113231);
	energy_match.check("total energy",        result.energy.total(),          -9.802353624482);
	energy_match.check("kinetic energy",      result.energy.kinetic(),         4.200403363300);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,       0.602423636472);
	energy_match.check("Hartree energy",      result.energy.hartree,           0.219187944078);
	energy_match.check("external energy",     result.energy.external,         -0.562807745696);
	energy_match.check("non-local energy",    result.energy.nonlocal,          1.427229177688);
	energy_match.check("XC energy",           result.energy.xc,               -4.767994250622);
	energy_match.check("XC density integral", result.energy.nvxc,             -4.900777046978);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,       0.0);
	
	return energy_match.fail();
	
}

