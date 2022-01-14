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
	
	utils::match energy_match(5.0e-6);

	std::vector<input::atom> geo;

	auto distance = 2.2; //a bit larger than experiment to check the force
	
	geo.push_back( "N" | math::vector3<double>(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | math::vector3<double>(0.0, 0.0,  0.5*distance));

	auto box = systems::box::orthorhombic(10.0_b, 10.0_b, 12.0_b).cutoff_energy(40.0_Ha);
	
	systems::ions ions(box, geo);

	systems::electrons electrons(env.dist(), ions, box, input::config{});
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::dft(),
																				input::scf::davidson() | input::scf::density_mixing() | input::scf::energy_tolerance(1e-7_Ha) | input::scf::calculate_forces());

	for(int iatom = 0; iatom < result.forces.size(); iatom++){
		printf("Force atom %d = %20.14f %20.14f %20.14f\n", iatom, result.forces[iatom][0], result.forces[iatom][1], result.forces[iatom][2]);
	}
	
	energy_match.check("ion-ion energy",      result.energy.ion,             -1.517434464849);

	energy_match.check("total energy",        result.energy.total(),         -20.642638450082);
	energy_match.check("kinetic energy",      result.energy.kinetic(),        13.163479187955);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,      -5.266545183839);
	energy_match.check("Hartree energy",      result.energy.hartree,          14.494239370982);
	energy_match.check("external energy",     result.energy.external,        -39.444396029142);
	energy_match.check("non-local energy",    result.energy.nonlocal,         -1.620128453793);
	energy_match.check("XC energy",           result.energy.xc,               -5.718398061234);
	energy_match.check("XC density integral", result.energy.nvxc,             -6.353978630822);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,       0.0);
	
	energy_match.check("force 1 x",           result.forces[0][0],            -0.000000025966);
	energy_match.check("force 1 y",           result.forces[0][1],             0.000000017231);
	energy_match.check("force 1 z",           result.forces[0][2],             0.145730384341);
	energy_match.check("force 2 x",           result.forces[1][0],             0.000000035122);
	energy_match.check("force 2 y",           result.forces[1][1],            -0.000000025695);
	energy_match.check("force 2 z",           result.forces[1][2],            -0.145731033675);

	//		EHRENFEST HASN'T BEEN TESTED YET
	
	auto dt = 0.025_atomictime;
	
	auto propagation = real_time::propagate(ions, electrons, input::interaction::dft(), input::rt::num_steps(10) | input::rt::dt(dt), ions::propagator::molecular_dynamics{});
		
	auto ofs = std::ofstream{"td_coordinates.dat"};
	
	for(unsigned long ii = 0; ii < propagation.coordinates.size(); ii++){
		ofs << propagation.time[ii]
				<< '\t' << propagation.coordinates[ii][0]
				<< '\t' << propagation.velocities[ii][0]
				<< '\t' << propagation.forces[ii][0] << std::endl;		
	}

	return energy_match.fail();

}

