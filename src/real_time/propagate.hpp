/* -*- indent-tabs-mode: t -*- */
#ifndef INQ__REAL_TIME__PROPAGATE
#define INQ__REAL_TIME__PROPAGATE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <hamiltonian/self_consistency.hpp>
#include <hamiltonian/forces.hpp>
#include <operations/overlap_diagonal.hpp>
#include <observables/dipole.hpp>
#include <options/real_time.hpp>
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
void propagate(systems::ions & ions, systems::electrons & electrons, ProcessFunction func, const options::theory & inter, const options::real_time & opts, IonSubPropagator const& ion_propagator = {}, Perturbation const & pert = {}){
		CALI_CXX_MARK_FUNCTION;
		
		auto console = electrons.logger();
		if(console) console->trace("initializing real-time propagation");

		const double dt = opts.dt();
		const int numsteps = opts.num_steps();

		for(auto & phi : electrons.kpin()) pert.zero_step(phi);
		
		electrons.spin_density() = observables::density::calculate(electrons);

		hamiltonian::self_consistency sc(inter, electrons.states_basis(), electrons.density_basis(), electrons.states().num_density_components(), pert);
		hamiltonian::ks_hamiltonian<complex> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(),
																						 ions, sc.exx_coefficient(), /* use_ace = */ opts.propagator() == options::real_time::electron_propagator::CRANK_NICOLSON);
		hamiltonian::energy energy;

		sc.update_ionic_fields(electrons.states_comm(), ions, electrons.atomic_pot());
		sc.update_hamiltonian(ham, energy, electrons.spin_density(), /* time = */ 0.0);

		ham.exchange.update(electrons);

		energy.calculate(ham, electrons);
		energy.ion(inq::ions::interaction_energy(ions.cell(), ions, electrons.atomic_pot()));

		auto forces = decltype(hamiltonian::calculate_forces(ions, electrons, ham)){};
		
		if(ion_propagator.needs_force) forces = hamiltonian::calculate_forces(ions, electrons, ham);

		func(real_time::viewables{false, 0, 0.0, ions, electrons, energy, forces, ham, pert});

		if(console) console->trace("starting real-time propagation");
		if(console) console->info("step {:9d} :  t =  {:9.3f}  e = {:.12f}", 0, 0.0, energy.total());

		auto iter_start_time = std::chrono::high_resolution_clock::now();
		for(int istep = 0; istep < numsteps; istep++){
			CALI_CXX_MARK_SCOPE("time_step");

			switch(opts.propagator()){
			case options::real_time::electron_propagator::ETRS :
				etrs(istep*dt, dt, ions, electrons, ion_propagator, forces, ham, sc, energy);
				break;
			case options::real_time::electron_propagator::CRANK_NICOLSON :
				crank_nicolson(istep*dt, dt, ions, electrons, ion_propagator, forces, ham, sc, energy);
				break;
			}

			energy.calculate(ham, electrons);
			
			if(ion_propagator.needs_force) forces = hamiltonian::calculate_forces(ions, electrons, ham);
			
			//propagate ionic velocities to t + dt
			ion_propagator.propagate_velocities(dt, ions, forces);

			func(real_time::viewables{istep == numsteps - 1, istep, (istep + 1.0)*dt, ions, electrons, energy, forces, ham, pert});			
			
			auto new_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = new_time - iter_start_time;
			
			if(console) console->info("step {:9d} :  t =  {:9.3f}  e = {:.12f}  wtime = {:9.3f}", istep + 1, (istep + 1)*dt, energy.total(), elapsed_seconds.count());

			iter_start_time = new_time;
		}

		if(console) console->trace("real-time propagation ended normally");
	}
}
}
#endif

#ifdef INQ_REAL_TIME_PROPAGATE_UNIT_TEST
#undef INQ_REAL_TIME_PROPAGATE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
