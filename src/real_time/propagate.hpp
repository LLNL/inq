/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2020-2021 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__REAL_TIME__PROPAGATE
#define INQ__REAL_TIME__PROPAGATE

#include <systems/ions.hpp>
#include <hamiltonian/calculate_energy.hpp>
#include <hamiltonian/self_consistency.hpp>
#include <hamiltonian/forces.hpp>
#include <operations/overlap_diagonal.hpp>
#include <observables/dipole.hpp>
#include <ions/propagator.hpp>
#include <systems/electrons.hpp>
#include <real_time/crank_nicolson.hpp>
#include <real_time/etrs.hpp>
#include <real_time/result.hpp>
#include <utils/profiling.hpp>

#include <chrono>

namespace inq {
namespace real_time {

template<typename IonSubPropagator = ions::propagator::fixed>
real_time::result propagate(systems::ions & ions, systems::electrons & electrons, const input::interaction & inter, const input::rt & options, IonSubPropagator const& ion_propagator = {}){

		CALI_CXX_MARK_FUNCTION;
		
		const double dt = options.dt();
		const int numsteps = options.num_steps();

		result res;

		electrons.density_ = density::calculate(electrons);

		hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_);
		hamiltonian::ks_hamiltonian<basis::real_space> ham(electrons.states_basis_, ions.cell(), electrons.atomic_pot_, inter.fourier_pseudo_value(), ions.geo(),
																											 electrons.states_.num_states(), sc.exx_coefficient(), electrons.states_basis_comm_);
		hamiltonian::energy energy;
		
		sc.update_ionic_fields(electrons.states_comm_, ions, electrons.atomic_pot_);
		
		ham.scalar_potential = sc.ks_potential(electrons.density_, energy);

		auto ecalc = hamiltonian::calculate_energy(ham, electrons);
		energy.eigenvalues = ecalc.sum_eigenvalues_;
		
		energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);
		
		if(electrons.full_comm_.root()) tfm::format(std::cout, "step %9d :  t =  %9.3f  e = %.12f\n", 0, 0.0, energy.total());

		auto forces = decltype(hamiltonian::calculate_forces(ions, electrons, ham)){};
		
		if(ion_propagator.needs_force) forces = hamiltonian::calculate_forces(ions, electrons, ham);

		res.save_iteration_results(0.0, ions, electrons, energy, forces);

		auto iter_start_time = std::chrono::high_resolution_clock::now();
		for(int istep = 0; istep < numsteps; istep++){
			CALI_CXX_MARK_SCOPE("time_step");

			//propagate using the chosen method
			switch(options.propagator()){
			case input::rt::electron_propagator::ETRS :
				etrs(dt, ions, electrons, ion_propagator, forces, ham, sc, energy);
				break;
			case input::rt::electron_propagator::CRANK_NICOLSON :
				crank_nicolson(dt, ions, electrons, ion_propagator, forces, ham, sc, energy);
				break;
			}
			
			//calculate the new density, energy, forces
			electrons.density_ = density::calculate(electrons);
			ham.scalar_potential = sc.ks_potential(electrons.density_, energy);

			auto ecalc = hamiltonian::calculate_energy(ham, electrons);
			energy.eigenvalues = ecalc.sum_eigenvalues_;
			
			if(ion_propagator.needs_force) forces = hamiltonian::calculate_forces(ions, electrons, ham);
			
			//propagate ionic velocities to t + dt
			ion_propagator.propagate_velocities(dt, ions, forces);

			res.save_iteration_results((istep + 1.0)*dt, ions, electrons, energy, forces);

			auto new_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = new_time - iter_start_time;
			if(electrons.full_comm_.root()) tfm::format(std::cout, "step %9d :  t =  %9.3f  e = %.12f  wtime = %9.3f\n", istep + 1, (istep + 1)*dt, energy.total(), elapsed_seconds.count());

			iter_start_time = new_time;
		}

		return res;
	}
}
}

#ifdef INQ_REAL_TIME_PROPAGATE_UNIT_TEST
#undef INQ_REAL_TIME_PROPAGATE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("real_time::propagate", "[real_time::propagate]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
