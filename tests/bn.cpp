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

  input::poscar vasp_file(config::path::unit_tests_data() + "bn.poscar");

	auto box = systems::box(vasp_file.lattice());
  
	systems::ions ions(box);
	
	ions.insert(vasp_file.atoms());
  
  systems::electrons electrons(env.par(), ions, box, input::config::cutoff(35.0_Ha) | input::config::extra_states(3), input::kpoints::grid({2, 2, 2}, true));
	
  ground_state::initial_guess(ions, electrons);
	
  auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
	
  energy_match.check("total energy",        result.energy.total(),         -13.415882371703);
  energy_match.check("kinetic energy",      result.energy.kinetic(),         9.561946750259);
  energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.946744603927);
  energy_match.check("Hartree energy",      result.energy.hartree(),         1.768191600138);
  energy_match.check("external energy",     result.energy.external(),       -7.925687577608);
  energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.120110167160);
  energy_match.check("XC energy",           result.energy.xc(),             -4.410160558293);
  energy_match.check("XC density integral", result.energy.nvxc(),           -4.999276809695);
  energy_match.check("ion-ion energy",      result.energy.ion(),           -11.290062419039);
	
	return energy_match.fail();
	
}


