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

	auto dcc = 1.42_A;
  auto aa = sqrt(3)*dcc;
  auto lz = 10.0_b;

	systems::ions ions(systems::cell::lattice(aa*vector3{1.0, 0.0, 0.0}, aa*vector3{-1.0/2.0, sqrt(3.0)/2.0, 0.0}, {0.0_b, 0.0_b, lz}).periodicity(2));
	
	ions.insert("C", {0.0_b, 0.0_b, 0.0_b});
	ions.insert("C", {0.0_b, dcc,   0.0_b});

	{
		systems::electrons electrons(env.par(), ions, options::electrons{}.spacing(aa/15.0).extra_states(2), input::kpoints::grid({1, 1, 1}, false));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.rpbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -11.805691503483);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         9.569622551441);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -3.508150621604);
		energy_match.check("Hartree energy",      result.energy.hartree(),       -11.111524294835);
		energy_match.check("external energy",     result.energy.external(),       15.082403332438);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.140029496383);
		energy_match.check("XC energy",           result.energy.xc(),             -4.392447131111);
		energy_match.check("XC density integral", result.energy.nvxc(),           -4.797098419430);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -19.813716465033);
		
	}
		
	fftw_cleanup();
	
	return energy_match.fail();
	
}


