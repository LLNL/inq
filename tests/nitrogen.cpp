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
	
	utils::match energy_match(6.0e-6);

	auto box = systems::box::orthorhombic(10.0_b, 10.0_b, 12.0_b).cutoff_energy(40.0_Ha);
	
	systems::ions ions(box);

	auto distance = 2.2_bohr; //a bit larger than experiment to check the force
	
	ions.insert("N", {0.0_b, 0.0_b, -0.5*distance});
	ions.insert("N", {0.0_b, 0.0_b,  0.5*distance});

	systems::electrons electrons(env.par(), ions, box, input::config{});
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::lda(), input::scf::energy_tolerance(1e-8_Ha) | input::scf::calculate_forces());

	for(int iatom = 0; iatom < result.forces.size(); iatom++){
		printf("Force atom %d = %20.14f %20.14f %20.14f\n", iatom, result.forces[iatom][0], result.forces[iatom][1], result.forces[iatom][2]);
	}
	
	energy_match.check("ion-ion energy",      result.energy.ion,             -1.517434464849);

	energy_match.check("total energy",        result.energy.total(),         -20.642638450082);
	energy_match.check("kinetic energy",      result.energy.kinetic(),        13.163479201031);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,      -5.266545173082);
	energy_match.check("Hartree energy",      result.energy.hartree,          14.494239381895);
	energy_match.check("external energy",     result.energy.external,        -39.444396051787);
	energy_match.check("non-local energy",    result.energy.nonlocal,         -1.620128452378);
	energy_match.check("XC energy",           result.energy.xc,               -5.718398063994);
	energy_match.check("XC density integral", result.energy.nvxc,             -6.353978633738);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,       0.0);
	
	energy_match.check("force 1 x",           result.forces[0][0],             0.00000000516220);
	energy_match.check("force 1 y",           result.forces[0][1],             0.00000001909934);
	energy_match.check("force 1 z",           result.forces[0][2],             0.14573044507129);
	energy_match.check("force 2 x",           result.forces[1][0],             0.00000001039003);
	energy_match.check("force 2 y",           result.forces[1][1],             0.00000002036321);
	energy_match.check("force 2 z",           result.forces[1][2],            -0.14573095235185);
	
	//		EHRENFEST HASN'T BEEN TESTED YET

	auto ofs = std::ofstream{"td_coordinates.dat"};
	auto process = [&ofs](auto data){
		ofs << data.time() << '\t' << data.coordinates(0) << '\t' << data.velocities(0) << '\t' << data.forces(0) << std::endl;		
	};
	
	auto dt = 0.025_atomictime;
	real_time::propagate(ions, electrons, process, input::interaction::lda(), input::rt::num_steps(10) | input::rt::dt(dt), ions::propagator::molecular_dynamics{});

	return energy_match.fail();

}

