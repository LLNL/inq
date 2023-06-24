/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);
	
	utils::match energy_match(6.0e-6);

	auto box = systems::box::orthorhombic(10.0_b, 10.0_b, 12.0_b);
	
	systems::ions ions(box);

	auto distance = 2.2_bohr; //a bit larger than experiment to check the force
	
	ions.insert("N", {0.0_b, 0.0_b, -0.5*distance});
	ions.insert("N", {0.0_b, 0.0_b,  0.5*distance});

	systems::electrons electrons(env.par(), ions, input::config::cutoff(40.0_Ha));
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::lda(), input::scf::energy_tolerance(1e-9_Ha) | input::scf::calculate_forces());

	for(int iatom = 0; iatom < result.forces.size(); iatom++){
		printf("Force atom %d = %20.14f %20.14f %20.14f\n", iatom, result.forces[iatom][0], result.forces[iatom][1], result.forces[iatom][2]);
	}
	
	energy_match.check("ion-ion energy",      result.energy.ion(),           -1.517434464849);

	energy_match.check("total energy",        result.energy.total(),         -20.642638450082);
	energy_match.check("kinetic energy",      result.energy.kinetic(),        13.163479200329);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -5.266545154146);
	energy_match.check("Hartree energy",      result.energy.hartree(),        14.494239402365);
	energy_match.check("external energy",     result.energy.external(),      -39.444396066271);
	energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.620128454760);
	energy_match.check("XC energy",           result.energy.xc(),             -5.718398066896);
	energy_match.check("XC density integral", result.energy.nvxc(),           -6.353978638174);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange(),     0.0);
	
	energy_match.check("force 1 x",           result.forces[0][0],            -0.00000000250910);
	energy_match.check("force 1 y",           result.forces[0][1],            -0.00000000185763);
	energy_match.check("force 1 z",           result.forces[0][2],             0.14573092924241);
	energy_match.check("force 2 x",           result.forces[1][0],             0.00000000341961);
	energy_match.check("force 2 y",           result.forces[1][1],             0.00000000029088);
	energy_match.check("force 2 z",           result.forces[1][2],            -0.14573050416125);
	
	//		EHRENFEST HASN'T BEEN TESTED YET

	std::ofstream ofs;
	if(electrons.root()) ofs.open("td_coordinates.dat");
	
	auto process = [&ofs](auto data){
		if(data.root()) ofs << data.time() << '\t' << data.coordinates(0) << '\t' << data.velocities(0) << '\t' << data.forces(0) << std::endl;		
	};
	
	auto dt = 0.025_atomictime;
	real_time::propagate(ions, electrons, process, input::interaction::lda(), input::rt::num_steps(10) | input::rt::dt(dt), ions::propagator::molecular_dynamics{});

	return energy_match.fail();

}

