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

#ifndef INQ__REAL_TIME__CRANK_NICOLSON
#define INQ__REAL_TIME__CRANK_NICOLSON

#include <density/calculate.hpp>
#include <operations/preconditioner.hpp>
#include <solvers/steepest_descent.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace real_time {

template <class IonSubPropagator, class ForcesType, class HamiltonianType, class SelfConsistencyType, class EnergyType>
void crank_nicolson(double const dt, systems::ions & ions, systems::electrons & electrons, IonSubPropagator const & ion_propagator, ForcesType const & forces, HamiltonianType & ham, SelfConsistencyType & sc, EnergyType & energy){

	auto op = [&](auto phi){
		auto opphi = ham(phi);

		gpu::run(opphi.fields().local_set_size(), opphi.fields().basis().local_size(),
						 [opph = begin(opphi.matrix()), ph = begin(phi.matrix()), dt] GPU_LAMBDA (auto ist, auto ip){
							 opph[ip][ist] = ph[ip][ist] + complex{0.0, 0.5*dt}*opph[ip][ist];
						 });

		return opphi;
	};


	auto prec = [&](auto phi){
	};

	std::vector<states::orbital_set<basis::real_space, complex>> rhs; 
	rhs.reserve(electrons.lot_size());
	
	//calculate the right hand side with H(t)
	auto iphi = 0;
	for(auto & phi : electrons.lot()){
		rhs.emplace_back(ham(phi));
		gpu::run(rhs[iphi].fields().local_set_size(), rhs[iphi].fields().basis().local_size(),
						 [rh = begin(rhs[iphi].matrix()), ph = begin(phi.matrix()), dt] GPU_LAMBDA (auto ist, auto ip){
							 rh[ip][ist] = (ph[ip][ist] - complex{0.0, 0.5*dt}*rh[ip][ist]);
						 });
		iphi++;
	}
	
	//propagate ionic positions to t + dt
	ion_propagator.propagate_positions(dt, ions, forces);
	if(not ion_propagator.static_ions) {
		sc.update_ionic_fields(ions, electrons.atomic_pot_);
		ham.update_projectors(electrons.states_basis_, ions.cell(), electrons.atomic_pot_, ions.geo());
		energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);
	}

	ham.scalar_potential = sc.ks_potential(electrons.density_, energy);

	//now calculate the wave functions in t + dt by solving a linear equation
	for(int istep = 0; istep < 10; istep++) {
		auto iphi = 0;
		for(auto & phi : electrons.lot()){
			solvers::steepest_descent(op, prec, rhs[iphi], phi);
			iphi++;
		}

	}
	
}

}
}

#ifdef INQ_REAL_TIME_CRANK_NICOLSON_UNIT_TEST
#undef INQ_REAL_TIME_CRANK_NICOLSON_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("real_time::crank_nicolson", "[real_time::crank_nicolson]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
