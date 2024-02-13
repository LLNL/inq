/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__PROPAGATOR
#define INQ__IONS__PROPAGATOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <solvers/velocity_verlet.hpp>

namespace inq {
namespace ions {
namespace propagator {

struct fixed {

	constexpr bool static_ions() const {
		return true;
	}

	constexpr bool needs_force() const {
		return false;
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons &, TypeForces const &){
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){
	}

};

struct impulsive {

	constexpr bool static_ions() const {
		return false;
	}

	constexpr bool needs_force() const {
		return false;
	}
	
	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons& ions, TypeForces const &){
		for(int i = 0; i != ions.size(); ++i)
			ions.positions()[i] += dt*ions.velocities()[i];
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){}

};


struct molecular_dynamics{

	constexpr bool static_ions() const {
		return false;
	}

	constexpr bool needs_force() const {
		return true;
	}

	template <typename TypeIons, typename TypeForces>
	static auto acceleration(TypeIons& ions, TypeForces forces){

		for(int iatom = 0; iatom < ions.size(); iatom++) forces[iatom] /= ions.atoms()[iatom].mass();
		return forces;

	}
	
	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons& ions, TypeForces const & forces){
		solvers::velocity_verlet::propagate_positions(dt, acceleration(ions, forces), ions.velocities(), ions.positions());
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons & ions, TypeForces const & forces){
		solvers::velocity_verlet::propagate_velocities(dt, acceleration(ions, forces), ions.velocities());
	}


};



}
}
}
#endif

#ifdef INQ_IONS_PROPAGATOR_UNIT_TEST
#undef INQ_IONS_PROPAGATOR_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
