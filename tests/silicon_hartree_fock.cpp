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
	
	ions.insert_fractional("Si", {0.0,  0.0,  0.0 });
	ions.insert_fractional("Si", {0.25, 0.25, 0.25});
	ions.insert_fractional("Si", {0.5,  0.5,  0.0 });
	ions.insert_fractional("Si", {0.75, 0.75, 0.25});
	ions.insert_fractional("Si", {0.5,  0.0,  0.5 });
	ions.insert_fractional("Si", {0.75, 0.25, 0.75});
	ions.insert_fractional("Si", {0.0,  0.5,  0.5 });
	ions.insert_fractional("Si", {0.25, 0.75, 0.75});

	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 4) parstates = 2;	
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	systems::electrons electrons(env.par().states(parstates), ions, box, input::config::extra_states(4), input::kpoints::grid({1, 1, 1}, true));

	ground_state::initial_guess(ions, electrons);
	ground_state::calculate(ions, electrons, input::interaction::dft(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-4_Ha));
 	auto result = ground_state::calculate(ions, electrons, input::interaction::hartree_fock(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
	
	energy_match.check("total energy",        result.energy.total(),       -30.495900907055);
	energy_match.check("kinetic energy",      result.energy.kinetic(),      13.234089410583);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),  -6.187488553120);
	energy_match.check("hartree energy",      result.energy.hartree(),       2.461774760792);
	energy_match.check("external energy",     result.energy.external(),     -9.170795269265);
	energy_match.check("non-local energy",    result.energy.nonlocal(),      4.099633587891);
	energy_match.check("XC energy",           result.energy.xc(),            0.000000000000);
	energy_match.check("XC density integral", result.energy.nvxc(),          0.000000000000);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange(),  -9.636982901957);
	energy_match.check("ion-ion energy",      result.energy.ion(),         -31.483620495100);
	
	electrons.save("silicon_restart");

	auto ked = observables::kinetic_energy_density(electrons);

	energy_match.check("kinetic energy", operations::integral(ked), 13.234089398748);
	
	fftw_cleanup(); //required for valgrind
	
	return energy_match.fail();
	
}

