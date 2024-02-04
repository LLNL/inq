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

	utils::match match(3.0e-5);

	systems::ions ions(systems::cell::cubic(10.0_b).finite());

	auto distance = 121.0_pm;
	
	ions.insert("O", {-distance/2, 0.0_b, 0.0_b});
	ions.insert("O", {distance/2, 0.0_b, 0.0_b});	

	{
		systems::electrons electrons(ions, options::electrons{}.spacing(0.43_b).spin_polarized().temperature(1000.0_K).extra_states(2));
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), options::ground_state{}.mixing(0.2).energy_tolerance(1e-8_Ha));
		
		match.check("total energy",        result.energy.total()    ,   -32.885878270495);
		match.check("kinetic energy",      result.energy.kinetic()   ,   20.663840649123);
		match.check("eigenvalues",         result.energy.eigenvalues(),  -7.271345986319);
		match.check("hartree energy",      result.energy.hartree(),      42.107413515539);
		match.check("external energy",     result.energy.external(),    -98.967748472806);
		match.check("non-local energy",    result.energy.nonlocal(),     -4.258411086657);
		match.check("XC energy",           result.energy.xc(),           -8.175088241373);
		match.check("XC density integral", result.energy.nvxc(),         -8.923854107057);
		match.check("HF exchange energy",  result.energy.exact_exchange(),   0.000000000000);
		match.check("ion-ion energy",      result.energy.ion(),          15.744115365679);
		
		auto magnetization = observables::magnetization(electrons.spin_density());
		auto total_mag = operations::integral(magnetization);
		
		match.check("magnetization x", total_mag[0], 0.0);
		match.check("magnetization y", total_mag[1], 0.0);
		match.check("magnetization z", total_mag[2], 2.0);

	}

	{
		systems::electrons electrons(ions, options::electrons{}.spacing(0.43_b).spin_non_collinear().temperature(1000.0_K).extra_states(2));
		ground_state::initial_guess(ions, electrons);
		ground_state::calculate(ions, electrons, options::theory{}.non_interacting(), options::ground_state{}.mixing(0.2).energy_tolerance(1e-8_Ha));
	}
		
	return match.fail();
}

