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
#include <perturbations/none.hpp>
#include <ions/propagator.hpp>
#include <systems/electrons.hpp>
#include <real_time/crank_nicolson.hpp>
#include <real_time/etrs.hpp>
#include <real_time/viewables.hpp>
#include <utils/profiling.hpp>

#include <chrono>

namespace inq {
namespace real_time {

template <typename ProcessFunction, typename IonSubPropagator = ions::propagator::fixed, typename Perturbation = perturbations::none>
void propagate(systems::ions & ions, systems::electrons & electrons, ProcessFunction func, const input::interaction & inter, const input::rt & options, IonSubPropagator const& ion_propagator = {}, Perturbation const & pert = {}){
		CALI_CXX_MARK_FUNCTION;
		
		const double dt = options.dt();
		const int numsteps = options.num_steps();

		for(auto & phi : electrons.lot()) pert.zero_step(phi);
		
		electrons.spin_density() = observables::density::calculate(electrons);

		hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_, electrons.states().num_density_components(), pert);
		hamiltonian::ks_hamiltonian<complex> ham(electrons.states_basis_, electrons.states(), electrons.atomic_pot_, inter.fourier_pseudo_value(), ions.geo(), electrons.states().num_states(), sc.exx_coefficient(), electrons.states_basis_comm_);
		hamiltonian::energy energy;

		sc.update_ionic_fields(electrons.states_comm_, ions, electrons.atomic_pot_);
		sc.update_hamiltonian(ham, energy, electrons.spin_density(), 0.0);

		auto ecalc = hamiltonian::calculate_energy(ham, electrons);
		energy.eigenvalues = ecalc.sum_eigenvalues_;
		
		energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);
		
		if(electrons.full_comm_.root()) tfm::format(std::cout, "step %9d :  t =  %9.3f  e = %.12f\n", 0, 0.0, energy.total());

		auto forces = decltype(hamiltonian::calculate_forces(ions, electrons, ham)){};
		
		if(ion_propagator.needs_force) forces = hamiltonian::calculate_forces(ions, electrons, ham);

		func(real_time::viewables{false, 0, 0.0, ions, electrons, energy, forces, ham, pert});
		
		auto iter_start_time = std::chrono::high_resolution_clock::now();
		for(int istep = 0; istep < numsteps; istep++){
			CALI_CXX_MARK_SCOPE("time_step");

			//propagate using the chosen method
			switch(options.propagator()){
			case input::rt::electron_propagator::ETRS :
				etrs(istep*dt, dt, ions, electrons, ion_propagator, forces, ham, sc, energy);
				break;
			case input::rt::electron_propagator::CRANK_NICOLSON :
				crank_nicolson(istep*dt, dt, ions, electrons, ion_propagator, forces, ham, sc, energy);
				break;
			}
			
			//calculate the new density, energy, forces
			electrons.spin_density() = observables::density::calculate(electrons);
			sc.update_hamiltonian(ham, energy, electrons.spin_density(), (istep + 1.0)*dt);
			
			auto ecalc = hamiltonian::calculate_energy(ham, electrons);
			energy.eigenvalues = ecalc.sum_eigenvalues_;
			
			if(ion_propagator.needs_force) forces = hamiltonian::calculate_forces(ions, electrons, ham);
			
			//propagate ionic velocities to t + dt
			ion_propagator.propagate_velocities(dt, ions, forces);

			func(real_time::viewables{istep == numsteps - 1, istep, (istep + 1.0)*dt, ions, electrons, energy, forces, ham, pert});			
			
			auto new_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = new_time - iter_start_time;
			if(electrons.full_comm_.root()) tfm::format(std::cout, "step %9d :  t =  %9.3f  e = %.12f  wtime = %9.3f\n", istep + 1, (istep + 1)*dt, energy.total(), elapsed_seconds.count());

			iter_start_time = new_time;
		}
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
