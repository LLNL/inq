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

#include <density/calculate.hpp>
#include <operations/exponential.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace real_time {

template <class IonSubPropagator, class ForcesType, class HamiltonianType, class SelfConsistencyType, class EnergyType>
void etrs(double const dt, systems::ions & ions, systems::electrons & electrons, IonSubPropagator const & ion_propagator, ForcesType const & forces, HamiltonianType & ham, SelfConsistencyType & sc, EnergyType & energy){
	
	electrons.density_ = 0.0;
	int iphi = 0;
	for(auto & phi : electrons.lot()){
		
		//propagate half step and full step with H(t)
		auto fullstep_phi = operations::exponential_2_for_1(ham, complex(0.0, dt), complex(0.0, dt/2.0), phi);
		
		//calculate H(t + dt) from the full step propagation
		density::calculate_add(electrons.occupations()[iphi], fullstep_phi, electrons.density_);

		iphi++;
	}

	if(electrons.lot_states_comm_.size() > 1){
		electrons.lot_states_comm_.all_reduce_in_place_n(raw_pointer_cast(electrons.density_.linear().data_elements()), electrons.density_.linear().size(), std::plus<>{});
	}

	//propagate ionic positions to t + dt
	ion_propagator.propagate_positions(dt, ions, forces);
	if(not ion_propagator.static_ions) {
		sc.update_ionic_fields(electrons.states_comm_, ions, electrons.atomic_pot_);
		ham.update_projectors(electrons.states_basis_, ions.cell(), electrons.atomic_pot_, ions.geo());
		energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);
	}

	ham.scalar_potential = sc.ks_potential(electrons.density_, energy);

	//propagate the other half step with H(t + dt)
	for(auto & phi : electrons.lot()){
		operations::exponential_in_place(ham, complex(0.0, dt/2.0), phi);
	}
	
}

}
}

#ifdef INQ_REAL_TIME_ETRS_UNIT_TEST
#undef INQ_REAL_TIME_ETRS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("real_time::etrs", "[real_time::etrs]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
