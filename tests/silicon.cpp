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

	input::environment env{};

	utils::match energy_match(3.0e-5);

	auto a = 10.18_b;
	systems::ions ions(systems::cell::cubic(a));
	
	ions.insert_fractional("Si", {0.0,  0.0,  0.0 });
	ions.insert_fractional("Si", {0.25, 0.25, 0.25});
	ions.insert_fractional("Si", {0.5,  0.5,  0.0 });
	ions.insert_fractional("Si", {0.75, 0.75, 0.25});
	ions.insert_fractional("Si", {0.5,  0.0,  0.5 });
	ions.insert_fractional("Si", {0.75, 0.25, 0.75});
	ions.insert_fractional("Si", {0.0,  0.5,  0.5 });
	ions.insert_fractional("Si", {0.25, 0.75, 0.75});

	systems::electrons electrons(env.par(), ions, input::kpoints::grid({2, 1, 1}, true), options::electrons{}.cutoff(25.0_Ha));

	ground_state::initial_guess(ions, electrons);

	auto result = ground_state::calculate(ions, electrons, options::theory{}.non_interacting(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-9_Ha));
	
	energy_match.check("total energy",     result.energy.total(),      -23.834202307770);
	energy_match.check("kinetic energy",   result.energy.kinetic(),     14.428062704088);
	energy_match.check("eigenvalues",      result.energy.eigenvalues(),  7.649418187330);
	energy_match.check("external energy",  result.energy.external(),   -12.019316573681);
	energy_match.check("non-local energy", result.energy.nonlocal(),     5.240672056923);
	energy_match.check("ion-ion energy",   result.energy.ion(),        -31.483620495100);
	
	auto ked = observables::kinetic_energy_density(electrons);

	energy_match.check("kinetic energy", operations::integral(ked), 14.428064504524);
	
	fftw_cleanup(); //required for valgrind
	
	return energy_match.fail();
	
}

