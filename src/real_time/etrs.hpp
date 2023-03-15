/* -*- indent-tabs-mode: t -*- */

/*
	Copyright (C) 2020-2022 Xavier Andrade, Alfredo A. Correa

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published by
	the Free Software Foundation; either version 3 of the License, or
	(at your option) any later version.
  
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.
  
	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef INQ__REAL_TIME__ETRS
#define INQ__REAL_TIME__ETRS

#include <observables/density.hpp>
#include <operations/exponential.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace real_time {

template <class IonSubPropagator, class ForcesType, class HamiltonianType, class SelfConsistencyType, class EnergyType>
void etrs(double const time, double const dt, systems::ions & ions, systems::electrons & electrons, IonSubPropagator const & ion_propagator, ForcesType const & forces, HamiltonianType & ham, SelfConsistencyType & sc, EnergyType & energy){

	int const nscf = 5;
	double const scf_threshold = 5e-5;

	systems::electrons::kpin_type save;
	
	int iphi = 0;
	for(auto & phi : electrons.kpin()){
		
		//propagate half step and full step with H(t)
		auto halfstep_phi = operations::exponential_2_for_1(ham, complex(0.0, dt/2.0), complex(0.0, dt), phi);
		save.emplace_back(std::move(halfstep_phi));
									 		
		iphi++;
	}

	electrons.spin_density() = observables::density::calculate(electrons);
	
	//propagate ionic positions to t + dt
	ion_propagator.propagate_positions(dt, ions, forces);
	if(not ion_propagator.static_ions) {
		sc.update_ionic_fields(electrons.states_comm(), ions, electrons.atomic_pot());
		ham.update_projectors(electrons.states_basis(), electrons.atomic_pot(), ions.geo());
		energy.ion(inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot()));
	}

	sc.update_hamiltonian(ham, energy, electrons.spin_density(), time + dt);
	ham.exchange.update(electrons);
	
	electrons.kpin() = save;
	
	//propagate the other half step with H(t + dt) self-consistently
	for(int iscf = 0; iscf < nscf; iscf++){

		int iphi = 0;
		for(auto & phi : electrons.kpin()) {
			if(iscf != 0) phi = save[iphi];
			operations::exponential_in_place(ham, complex(0.0, dt/2.0), phi);
			iphi++;
		}
		
		auto old_density = electrons.spin_density();
		electrons.spin_density() = observables::density::calculate(electrons);

		double delta = operations::integral_sum_absdiff(old_density, electrons.spin_density());
		auto done = (delta < scf_threshold) or (iscf == nscf - 1);
		
		sc.update_hamiltonian(ham, energy, electrons.spin_density(), time + dt);
		ham.exchange.update(electrons);
		
		if(done) break;
	}
	
}

}
}
#endif

#ifdef INQ_REAL_TIME_ETRS_UNIT_TEST
#undef INQ_REAL_TIME_ETRS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
