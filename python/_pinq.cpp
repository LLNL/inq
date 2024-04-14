/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#define INQ_PYTHON_INTERFACE
#include <inq/inq.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace inq;
using namespace inq::magnitude;

struct calculator {

private:
	
	options::theory theo_;
	options::electrons els_;
	ground_state::results results_;
	std::optional<systems::electrons> electrons_;
	
public:
	
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
		return results_.energy.total()*1.0_Ha/1.0_eV;
	}

	///////////////////////////////////
	
	auto get_forces(py::object atoms){

		py::array_t<double, py::array::c_style> forces_array({results_.forces.size(), 3l});
		
    auto arr = forces_array.mutable_unchecked();
		
    for (py::ssize_t iatom = 0; iatom < arr.shape(0); iatom++) {
			for (py::ssize_t idir = 0; idir < arr.shape(1); idir++) {
				arr(iatom, idir) = results_.forces[iatom][idir]*(1.0_Ha/1.0_eV)*(1.0_A/1.0_bohr); //convert to eV/A
			}
    }
		
    return forces_array;
	}

	///////////////////////////////////
	
	auto get_density(){

		assert(electrons_.has_value());

		auto const & dens = electrons_->spin_density().hypercubic();
			
		auto dx = std::get<0>(sizes(dens));
		auto dy = std::get<1>(sizes(dens));
		auto dz = std::get<2>(sizes(dens));
		auto dspin = std::get<3>(sizes(dens));		
		
		py::array_t<double, py::array::c_style> density_array({dx, dy, dz, dspin});
		
    auto arr = density_array.mutable_unchecked();
		
    for (py::ssize_t ix = 0; ix < arr.shape(0); ix++) {
			for (py::ssize_t iy = 0; iy < arr.shape(1); iy++) {
				for (py::ssize_t iz = 0; iz < arr.shape(2); iz++) {
					for (py::ssize_t ispin = 0; ispin < arr.shape(3); ispin++) {
						arr(ix, iy, iz, ispin) = dens[ix][iy][iz][ispin];
					}
				}
			}
    }
		
    return density_array;
	}

	///////////////////////////////////
	
	void calculate(py::object atoms){
		
		auto ions = systems::ions::import_ase(atoms);

		electrons_.emplace(systems::electrons(ions, els_));
		ground_state::initial_guess(ions, *electrons_);

		results_ = ground_state::calculate(ions, *electrons_, theo_, options::ground_state{}.energy_tolerance(1e-9_Ha).calculate_forces());
	}

	///////////////////////////////////
	
	auto scf_step(py::object atoms, py::array_t<double> const & potential){
		
		auto ions = systems::ions::import_ase(atoms);

		if(not electrons_.has_value()){
			electrons_.emplace(systems::electrons(ions, els_));
			ground_state::initial_guess(ions, *electrons_);
		}

		solvers::poisson poisson_solver;
		
		auto ionic_long_range = poisson_solver(electrons_->atomic_pot().ionic_density(electrons_->kpin_states_comm(), electrons_->density_basis(), ions));
		auto ionic_short_range = electrons_->atomic_pot().local_potential(electrons_->kpin_states_comm(), electrons_->density_basis(), ions);
		auto vion = operations::add(ionic_long_range, ionic_short_range);
		
		auto ham = hamiltonian::ks_hamiltonian<double>(electrons_->states_basis(), electrons_->brillouin_zone(), electrons_->states(), electrons_->atomic_pot(),
																									 ions, /*exx_coefficient = */ 0.0, /* use_ace = */ true);

    auto pot = potential.unchecked();
		
		for (py::ssize_t ix = 0; ix < pot.shape(0); ix++) {
			for (py::ssize_t iy = 0; iy < pot.shape(1); iy++) {
				for (py::ssize_t iz = 0; iz < pot.shape(2); iz++) {
					for (py::ssize_t ispin = 0; ispin < pot.shape(3); ispin++) {
						ham.scalar_potential().hypercubic()[ix][iy][iz][ispin] = vion.cubic()[ix][iy][iz] + pot(ix, iy, iz, ispin);
					}
				}
			}
    }

		//subspace diagonalization
		int ilot = 0;
		for(auto & phi : electrons_->kpin()) {
			electrons_->eigenvalues()[ilot] = ground_state::subspace_diagonalization(ham, phi);
			ilot++;
		}

		electrons_->update_occupations(electrons_->eigenvalues());

		operations::preconditioner prec;
		
		for(auto & phi : electrons_->kpin()) {
			auto fphi = operations::transform::to_fourier(std::move(phi));
			eigensolvers::steepest_descent(ham, prec, fphi);
			phi = operations::transform::to_real(std::move(fphi));
		}

		electrons_->spin_density() = observables::density::calculate(*electrons_);

		results_.energy.calculate(ham, *electrons_);
		
		return get_density();
	}
	
};

void clear() {
	interface::clear();
}

PYBIND11_MODULE(_pinq, module) {

	module.doc() = "Python interface for the INQ DFT/TDDFT library";

	py::class_<calculator>(module, "calculator")
		.def(py::init<py::args, py::kwargs&>())
		.def("get_potential_energy", &calculator::get_potential_energy)
		.def("get_forces",           &calculator::get_forces)
		.def("get_density",          &calculator::get_density)		
		.def("calculate",            &calculator::calculate)
		.def("scf_step",             &calculator::scf_step);

	auto interface_module = module.def_submodule("interface");
	interface::clear.python_interface(interface_module);
	interface::cell.python_interface(interface_module);
	interface::ions.python_interface(interface_module);
	
}
