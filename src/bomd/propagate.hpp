/* -*- indent-tabs-mode: t -*- */
#ifndef INQ__BOMD__PROPAGATE
#define INQ__BOMD__PROPAGATE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <ground_state/calculate.hpp>
#include <ions/propagator.hpp>
#include <options/real_time.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>
#include <utils/profiling.hpp>

#include <chrono>

namespace inq {
namespace bomd {

template <typename IonSubPropagator = ions::propagator::molecular_dynamics>
void propagate(systems::ions & ions, systems::electrons & electrons, const options::theory & inter, const options::real_time & opts, IonSubPropagator const& ion_propagator = {}){
		CALI_CXX_MARK_FUNCTION;
		
		auto console = electrons.logger();

		const double dt = opts.dt();
		const int numsteps = opts.num_steps();

		if(console) console->trace(std::string("initializing real-time propagation:\n") +
					   std::string("  time step        = {} atomictime ({:.2f} as)\n") + 
					   std::string("  number of steps  = {}\n") + 
					   std::string("  propagation time = {} atomictime ({:.2f} fs)"), dt, dt/0.041341373, numsteps, numsteps*dt, numsteps*dt/41.341373);


		auto res = ground_state::calculate(ions, electrons, inter, options::ground_state{}.calculate_forces());
		auto energy = res.energy;
		auto forces = res.forces;

		if(console) console->trace("starting Born-Oppenheimer propagation");
		if(console) console->info("step {:9d} :  t =  {:9.3f}  e = {:.12f}", 0, 0.0, energy.total());
		
		auto iter_start_time = std::chrono::high_resolution_clock::now();
		for(int istep = 0; istep < numsteps; istep++){
			CALI_CXX_MARK_SCOPE("time_step");

			//propagate ionic positions to t + dt
			ion_propagator.propagate_positions(dt, ions, forces);
			
			//calculate the electronic state at t + dt
			auto res = ground_state::calculate(ions, electrons, inter, options::ground_state{}.calculate_forces());
			energy = res.energy;
			forces = res.forces;

			//propagate ionic velocities to t + dt
			ion_propagator.propagate_velocities(dt, ions, forces);

			auto new_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = new_time - iter_start_time;
			
			if(console) console->info("step {:9d} :  t =  {:9.3f}  e = {:.12f}  wtime = {:9.3f}", istep + 1, (istep + 1)*dt, energy.total(), elapsed_seconds.count());

			iter_start_time = new_time;
		}

		if(console) console->trace("Born-Oppenheimer propagation ended normally");
	}

}
}
#endif

#ifdef INQ_BOMD_PROPAGATE_UNIT_TEST
#undef INQ_BOMD_PROPAGATE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
