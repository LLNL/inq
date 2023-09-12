/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;

auto ase_atoms_to_inq_ions(py::object atoms){

	using namespace inq;
	using namespace inq::magnitude;
	
	auto atomic_numbers = atoms.attr("get_atomic_numbers")().cast<py::array_t<int>>();
	auto num = static_cast<int *>(atomic_numbers.request().ptr);
	auto positions = atoms.attr("get_positions")().cast<py::array_t<double>>();
	auto pos = static_cast<double *>(positions.request().ptr);

	systems::ions ions(systems::cell::orthorhombic(10.0_b, 10.0_b, 12.0_b));

	for(int ii = 0; ii < atomic_numbers.size(); ii++){
		ions.insert(num[ii], 1.0_A*vector3{pos[3*ii + 0], pos[3*ii + 1], pos[3*ii + 2]});
	}

	return ions;
}

auto run(py::object atoms){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env{};
	
	utils::match energy_match(6.0e-6);

	auto ions = ase_atoms_to_inq_ions(atoms);
	
	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(40.0_Ha));
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), options::ground_state{}.energy_tolerance(1e-9_Ha).calculate_forces());

	return result.energy.total();
	
}

PYBIND11_MODULE(pinq, module) {

	module.doc() = "Python interface for the INQ DFT/TDDFT library";
	module.def("run", &run, "A function that runs a calculation");
	
}
