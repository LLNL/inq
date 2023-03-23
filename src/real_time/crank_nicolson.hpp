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

#include <observables/density.hpp>
#include <operations/preconditioner.hpp>
#include <solvers/steepest_descent.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace real_time {

template <class OpType>
struct crank_nicolson_op {
	OpType & op;
	complex factor;
	
	template <class PhiType>
	PhiType operator()(PhiType const & phi) const {
		auto opphi = op(phi);
		
		gpu::run(opphi.local_set_size(), opphi.basis().local_size(),
						 [opph = begin(opphi.matrix()), ph = begin(phi.matrix()), fac = factor] GPU_LAMBDA (auto ist, auto ip){
							 opph[ip][ist] = ph[ip][ist] + fac*opph[ip][ist];
						 });
		
		return opphi;
	}
};

template <class IonSubPropagator, class ForcesType, class HamiltonianType, class SelfConsistencyType, class EnergyType>
void crank_nicolson(double const time, double const dt, systems::ions & ions, systems::electrons & electrons, IonSubPropagator const & ion_propagator, ForcesType const & forces, HamiltonianType & ham, SelfConsistencyType & sc, EnergyType & energy){

	crank_nicolson_op<decltype(ham)> op{ham, complex{0.0, 0.5*dt}};
	crank_nicolson_op<decltype(ham)> op_rhs{ham, complex{0.0, -0.5*dt}};

	auto const st_tol = 1e-10;
	auto const dens_tol = 1e-5;	

	//calculate the right hand side with H(t)
	std::vector<states::orbital_set<basis::real_space, complex>> rhs; 
	rhs.reserve(electrons.kpin_size());	
	for(auto & phi : electrons.kpin()) rhs.emplace_back(op_rhs(phi));
	
	//propagate ionic positions to t + dt
	ion_propagator.propagate_positions(dt, ions, forces);
	if(not ion_propagator.static_ions) {
		sc.update_ionic_fields(electrons.states_comm(), ions, electrons.atomic_pot());
		ham.update_projectors(electrons.states_basis(), electrons.atomic_pot(), ions.geo());
		energy.ion(inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot()));
	}

	using mix_arr_type = std::remove_reference_t<decltype(electrons.spin_density().matrix().flatted())>;
	auto mixer = mixers::broyden<mix_arr_type>(4, 0.3, electrons.spin_density().matrix().flatted().size(), electrons.density_basis().comm());
	
	//now calculate the wave functions in t + dt by solving a self-consistent linear equation
	for(int istep = 0; istep < 200; istep++) {
		CALI_CXX_MARK_SCOPE("crank_nicolson:iteration");

		sc.update_hamiltonian(ham, energy, electrons.spin_density(), time + dt);
		ham.exchange.update(electrons);

		auto res = 0.0;
		auto iphi = 0;
		for(auto & phi : electrons.kpin()){
			auto ires = solvers::steepest_descent(op, operations::no_preconditioner{}, rhs[iphi], phi);
			res += ires*electrons.kpin_weights()[iphi];
			iphi++;
		}

		if(electrons.kpin_states_comm().size() > 1) electrons.kpin_states_comm().all_reduce_n(&res, 1, std::plus<>{});
		
		auto new_density = observables::density::calculate(electrons);
		auto density_diff = operations::integral_sum_absdiff(electrons.spin_density(), new_density)/electrons.states().num_electrons();
		
		std::cout << istep << '\t' << density_diff << '\t' << res << std::endl;
		
		if(res < st_tol and density_diff < dens_tol) break;
		
		auto tmp = +electrons.spin_density().matrix().flatted();
		mixer(tmp, new_density.matrix().flatted());
		electrons.spin_density().matrix().flatted() = tmp;
		observables::density::normalize(electrons.spin_density(), electrons.states().num_electrons());
	}

	electrons.spin_density() = observables::density::calculate(electrons);
	sc.update_hamiltonian(ham, energy, electrons.spin_density(), time + dt);
	ham.exchange.update(electrons);
	 
}

}
}
#endif

#ifdef INQ_REAL_TIME_CRANK_NICOLSON_UNIT_TEST
#undef INQ_REAL_TIME_CRANK_NICOLSON_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
