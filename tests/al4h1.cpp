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

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {1, boost::mpi3::fill});
	
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
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), input::scf::calculate_forces());

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

	energy_match.check("force 1 x",           result.forces[0][0],            -0.010958823420);
	energy_match.check("force 1 y",           result.forces[0][1],            -0.020170976619);
	energy_match.check("force 1 z",           result.forces[0][2],            -0.025825339342);
	energy_match.check("force 2 x",           result.forces[1][0],            -0.006997173079);
	energy_match.check("force 2 y",           result.forces[1][1],             0.028861371259);
	energy_match.check("force 2 z",           result.forces[1][2],             0.022386977302);
	energy_match.check("force 3 x",           result.forces[2][0],             0.021742257288);
	energy_match.check("force 3 y",           result.forces[2][1],            -0.015045246329);
	energy_match.check("force 3 z",           result.forces[2][2],             0.017083833900);
	energy_match.check("force 4 x",           result.forces[3][0],             0.018522478479);
	energy_match.check("force 4 y",           result.forces[3][1],             0.015573712792);
	energy_match.check("force 4 z",           result.forces[3][2],            -0.013419282975);
	energy_match.check("force 5 x",           result.forces[4][0],             0.032491963346);
	energy_match.check("force 5 y",           result.forces[4][1],            -0.014763643737);
	energy_match.check("force 5 z",           result.forces[4][2],             0.014730728104);
	
	return energy_match.fail();
	
}

