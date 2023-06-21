/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::input;
	using namespace inq::magnitude;
	
	environment env(argc, argv);

	auto distance = 1.06_angstrom;

	auto box = systems::box::orthorhombic(10.0_b, 10.0_b, 12.0_b).finite();
	
	systems::ions ions(box);

	ions.insert("N", {0.0_b, 0.0_b, -distance/2});
  ions.insert("N", {0.0_b, 0.0_b,  distance/2});
	
	systems::electrons electrons(env.par(), ions, box, input::config::cutoff(40.0_Ha));
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, interaction::pbe());

	std::cout << "N2 energy = " << result.energy.total() << std::endl;

}

