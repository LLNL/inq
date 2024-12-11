/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq::magnitude;

	inq::utils::match energy_match(3.0e-5);

	inq::systems::ions ions(inq::systems::cell::cubic(10.0_b).finite());
	ions.species_list().insert(inq::ionic::species("He").symbol("He_fr").pseudo_file(inq::config::path::unit_tests_data() + "He_fr.upf.gz"));
	ions.insert("He_fr", {0.0_b, 0.0_b, 0.0_b});
	
	inq::systems::electrons electrons(ions, inq::options::electrons{}.cutoff(40.0_Ha).spin_non_collinear().extra_states(2));
	inq::ground_state::initial_guess(ions, electrons);
	
	auto result = inq::ground_state::calculate(ions, electrons, inq::options::theory{}.pbe(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));

	/*
	energy_match.check("total energy",        result.energy.total(),          -0.445160072256);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         0.414315464604);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.234029035766);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.281309132025);
		energy_match.check("external energy",     result.energy.external(),       -0.909266382097);
		energy_match.check("non-local energy",    result.energy.non_local(),       0.000000000000);
		energy_match.check("XC energy",           result.energy.xc(),             -0.231518286787);
		energy_match.check("XC density integral", result.energy.nvxc(),           -0.301696382322);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),             0.000000000000);
	*/
	
	return energy_match.fail();
	
}
