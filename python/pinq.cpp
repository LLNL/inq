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

using namespace inq;
using namespace inq::magnitude;

auto ase_atoms_to_inq_ions(py::object atoms){

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

	options::theory theo_;
	options::electrons els_;
	ground_state::result result_;

	calculator(py::args const &, py::kwargs const & kwargs){
		auto const args_map = py::cast<std::unordered_map<std::string, py::object>>(kwargs);

		if(args_map.find("xc") != args_map.end()){
			auto functional = py::cast<std::string>(args_map.at("xc"));

			if(functional == "LDA" or functional == "lda") {
				theo_ = theo_.lda();
			} else if(functional == "PBE" or functional == "pbe") {
				theo_ = theo_.pbe();
			} else {
				throw std::runtime_error("pinq: Unknown functional '" + functional + "'.");
			}
		}

		if(args_map.find("ecut") != args_map.end()){
			els_ = els_.cutoff(1.0_Ry*py::cast<double>(args_map.at("ecut")));
		} else {
			throw std::runtime_error("pinq: Missing argument 'ecut'.");
		}

		if(args_map.find("extra_bands") != args_map.end()){
			els_ = els_.extra_states(py::cast<int>(args_map.at("extra_bands")));
		}
		
		if(args_map.find("width") != args_map.end()){
			els_ = els_.temperature(1.0_eV*py::cast<double>(args_map.at("width")));
		}
		
	}

	///////////////////////////////////
	
	auto get_potential_energy(py::object atoms){
		return result_.energy.total()*1.0_Ha/1.0_eV;
	}

	///////////////////////////////////
	
	auto get_forces(py::object atoms){

		py::array_t<double, py::array::c_style> forces_array({result_.forces.size(), 3l});
		
    auto arr = forces_array.mutable_unchecked();
		
    for (py::ssize_t iatom = 0; iatom < arr.shape(0); iatom++) {
			for (py::ssize_t idir = 0; idir < arr.shape(1); idir++) {
				arr(iatom, idir) = result_.forces[iatom][idir]*(1.0_Ha/1.0_eV)*(1.0_A/1.0_bohr); //convert to eV/A
			}
    }
		
    return forces_array;
	}

	///////////////////////////////////
	
	void calculate(py::object atoms){
		input::environment env{};
		
		utils::match energy_match(6.0e-6);
		
		auto ions = ase_atoms_to_inq_ions(atoms);

		systems::electrons electrons(env.par(), ions, els_);
		ground_state::initial_guess(ions, electrons);

		result_ = ground_state::calculate(ions, electrons, theo_, options::ground_state{}.energy_tolerance(1e-9_Ha).calculate_forces());

		std::cout << result_.forces[0] << std::endl;
		std::cout << result_.forces[1] << std::endl;		

	}
	
};

PYBIND11_MODULE(pinq, module) {

	module.doc() = "Python interface for the INQ DFT/TDDFT library";

	py::class_<calculator>(module, "calculator")
		.def(py::init<py::args, py::kwargs&>())
		.def("get_potential_energy", &calculator::get_potential_energy)
		.def("get_forces", &calculator::get_forces)
		.def("calculate", &calculator::calculate);

	
}
