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
	
	auto lattice = atoms.attr("get_cell")().attr("__array__")().cast<py::array_t<double>>();
	auto lat = static_cast<double *>(lattice.request().ptr);

	auto lat0 = vector3(1.0_A*lat[0], 1.0_A*lat[1], 1.0_A*lat[2]);
	auto lat1 = vector3(1.0_A*lat[3], 1.0_A*lat[4], 1.0_A*lat[5]);
	auto lat2 = vector3(1.0_A*lat[6], 1.0_A*lat[7], 1.0_A*lat[8]);
	
	systems::ions ions(systems::cell::lattice(lat0, lat1, lat2));

	auto atomic_numbers = atoms.attr("get_atomic_numbers")().cast<py::array_t<int>>();
	auto num = static_cast<int *>(atomic_numbers.request().ptr);
	auto positions = atoms.attr("get_positions")().cast<py::array_t<double>>();
	auto pos = static_cast<double *>(positions.request().ptr);

	for(int ii = 0; ii < atomic_numbers.size(); ii++){
		ions.insert(num[ii], 1.0_A*vector3{pos[3*ii + 0], pos[3*ii + 1], pos[3*ii + 2]});
	}

	return ions;
}

struct calculator {
	
	auto get_potential_energy(py::object atoms){
		
		using namespace inq;
		using namespace inq::magnitude;
		
		input::environment env{};
		
		utils::match energy_match(6.0e-6);
		
		auto ions = ase_atoms_to_inq_ions(atoms);
		
		systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(35.0_Ha));
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), options::ground_state{}.energy_tolerance(1e-9_Ha));
		
		return result.energy.total();
	}
	
};

PYBIND11_MODULE(pinq, module) {

	module.doc() = "Python interface for the INQ DFT/TDDFT library";

	py::class_<calculator>(module, "calculator")
		.def(py::init<>())
		.def("get_potential_energy", &calculator::get_potential_energy);
	
}
