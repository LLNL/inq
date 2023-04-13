/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <observables/kinetic_energy_density.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);

	utils::match energy_match(3.0e-5);

	std::vector<input::atom> geo;

	auto a = 10.18_b;

	auto box = systems::box::cubic(a).cutoff_energy(25.0_Ha);
	
	systems::ions ions(box);
	
	ions.insert("Si", {0.0_crys,  0.0_crys,  0.0_crys });
	ions.insert("Si", {0.25_crys, 0.25_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.5_crys,  0.0_crys });
	ions.insert("Si", {0.75_crys, 0.75_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.0_crys,  0.5_crys });
	ions.insert("Si", {0.75_crys, 0.25_crys, 0.75_crys});
	ions.insert("Si", {0.0_crys,  0.5_crys,  0.5_crys });
	ions.insert("Si", {0.25_crys, 0.75_crys, 0.75_crys});
	
	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	systems::electrons electrons(env.par().states(parstates), ions, box, input::config::extra_states(4), input::kpoints::grid({1, 1, 1}, true));

	ground_state::initial_guess(ions, electrons);
	ground_state::calculate(ions, electrons, input::interaction::dft(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-4_Ha));
 	auto result = ground_state::calculate(ions, electrons, input::interaction::hartree_fock(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
	
	energy_match.check("total energy",     result.energy.total()      , -30.495900918097);
	energy_match.check("kinetic energy",   result.energy.kinetic()    ,  13.234089373454);
	energy_match.check("eigenvalues",      result.energy.eigenvalues(),  -6.187488598514);
	energy_match.check("external energy",  result.energy.external()   ,  -9.170795172150);
	energy_match.check("non-local energy", result.energy.nonlocal()   ,   4.099633551216);
	energy_match.check("hf exchange",      result.energy.hf_exchange(),  -9.636982906989);
	energy_match.check("ion-ion energy",   result.energy.ion()        , -31.483620495100);
	
	electrons.save("silicon_restart");

	auto ked = observables::kinetic_energy_density(electrons);

	energy_match.check("kinetic energy", operations::integral(ked), 13.234089398748);
	
	fftw_cleanup(); //required for valgrid
	
	return energy_match.fail();
	
}

