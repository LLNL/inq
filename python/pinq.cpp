/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void run(){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env{};
	
	utils::match energy_match(6.0e-6);

	systems::ions ions(systems::cell::orthorhombic(10.0_b, 10.0_b, 12.0_b));

	auto distance = 2.2_bohr; //a bit larger than experiment to check the force
	
	ions.insert("N", {0.0_b, 0.0_b, -0.5*distance});
	ions.insert("N", {0.0_b, 0.0_b,  0.5*distance});

	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(40.0_Ha));
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), options::ground_state{}.energy_tolerance(1e-9_Ha).calculate_forces());

	std::cout << "Energy = " << result.energy.total() << std::endl;
	
}

PYBIND11_MODULE(pinq, module) {

	module.doc() = "Python interface for the INQ DFT/TDDFT library";
	module.def("run", &run, "A function that runs a calculation");
	
}
