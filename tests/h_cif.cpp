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
	
	auto ions = systems::ions::parse(config::path::unit_tests_data() + "H.cif");

	systems::electrons electrons(ions, options::electrons{}.cutoff(30.0_Ha));
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));

	utils::match energy_match(3.0e-5);

  energy_match.check("total energy",        result.energy.total(),          -2.355128752518);
  energy_match.check("kinetic energy",      result.energy.kinetic(),         2.125707369840);
  energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -1.511390321649);
  energy_match.check("Hartree energy",      result.energy.hartree(),         0.903321603071);
  energy_match.check("external energy",     result.energy.external(),       -3.319351558502);
  energy_match.check("non-local energy",    result.energy.non_local(),       -0.389801056728);
  energy_match.check("XC energy",           result.energy.xc(),             -1.351311210649);
  energy_match.check("XC density integral", result.energy.nvxc(),           -1.734588282400);
  energy_match.check("ion-ion energy",      result.energy.ion(),            -0.323693899550);

	auto all_eigenvalues = parallel::gather(+electrons.eigenvalues().flatted(), electrons.kpin_states_part(), electrons.kpin_states_comm(), 0);
	
	if(electrons.kpin_states_comm().root()){
		energy_match.check("eigenvalue 1",          all_eigenvalues[0],  -0.401182208666);
		energy_match.check("eigenvalue 2",          all_eigenvalues[1],  -0.354512952159);
	}
	
}
