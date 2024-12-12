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
	ions.species_list().pseudopotentials() = pseudo::set_id::pseudodojo_rel_pbe();
	ions.insert("He", {0.0_b, 0.0_b, 0.0_b});
	
	inq::systems::electrons electrons(ions, inq::options::electrons{}.cutoff(40.0_Ha).spin_non_collinear().extra_states(2));
	inq::ground_state::initial_guess(ions, electrons);
	
	auto result = inq::ground_state::calculate(ions, electrons, inq::options::theory{}.pbe(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));

	energy_match.check("total energy",        result.energy.total(),          -2.860762587218);
	energy_match.check("kinetic energy",      result.energy.kinetic(),         2.504793796156);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -1.161792830932);
	energy_match.check("Hartree energy",      result.energy.hartree(),         1.986218327002);
	energy_match.check("external energy",     result.energy.external(),       -4.892155878650);
	energy_match.check("non-local energy",    result.energy.non_local(),      -1.440150012679);
	energy_match.check("XC energy",           result.energy.xc(),             -1.019468819047);
	energy_match.check("XC density integral", result.energy.nvxc(),           -1.306717389763);
	energy_match.check("HF exchange energy",  result.energy.exact_exchange(),  0.000000000000);
	energy_match.check("ion-ion energy",      result.energy.ion(),             0.000000000000);
	
	return energy_match.fail();
	
}
