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

	utils::match energy_match(3.0e-5);

  auto ions = systems::ions::parse(config::path::unit_tests_data() + "bn.poscar");

  systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(35.0_Ha).extra_states(3), input::kpoints::grid({2, 2, 2}, true));
	
  ground_state::initial_guess(ions, electrons);
	
  auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
	
  energy_match.check("total energy",        result.energy.total(),         -13.359343074445);
  energy_match.check("kinetic energy",      result.energy.kinetic(),         9.324415990405);
  energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.921653047910);
  energy_match.check("Hartree energy",      result.energy.hartree(),         1.720766003426);
  energy_match.check("external energy",     result.energy.external(),       -7.771139056402);
  energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.950079325629);
  energy_match.check("XC energy",           result.energy.xc(),             -4.393244267206);
  energy_match.check("XC density integral", result.energy.nvxc(),           -4.966382663135);
  energy_match.check("ion-ion energy",      result.energy.ion(),           -11.290062419039);
	
	return energy_match.fail();
	
}


