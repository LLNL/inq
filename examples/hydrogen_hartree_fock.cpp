/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	inq::input::environment env(argc, argv);
	
	inq::utils::match energy_match(3.0e-5);

	auto local_h = inq::input::species("H").symbol("Hloc").pseudo(inq::config::path::unit_tests_data() + "H.blyp-vbc.UPF");
	
	inq::systems::ions ions(inq::systems::cell::cubic(15.0_b).finite());
	ions.insert(local_h, {150.0_b, -30.0_b, 0.0_b});

	inq::systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(40.0_Ha));
	inq::ground_state::initial_guess(ions, electrons);
	
	inq::ground_state::calculate(ions, electrons);
	auto result = inq::ground_state::calculate(ions, electrons, inq::options::theory{}.hartree_fock(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));
	
	energy_match.check("total energy",        result.energy.total(),      -0.578525486338);
	energy_match.check("kinetic energy",      result.energy.kinetic(),     0.348185715818);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),-0.230237311450);
	energy_match.check("Hartree energy",      result.energy.hartree(),     0.254438812760);
	energy_match.check("external energy",     result.energy.external(),   -0.832861830904);
	energy_match.check("non-local energy",    result.energy.nonlocal(),    0.0);
	energy_match.check("XC energy",           result.energy.xc(),          0.0);
	energy_match.check("XC density integral", result.energy.nvxc(),        0.0);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange(),-0.254438821884);
	energy_match.check("ion-ion energy",      result.energy.ion(),        -0.093849362128);
	
	return energy_match.fail();
	
}
