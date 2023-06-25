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

	static constexpr bool static_ions = true;
	static constexpr bool needs_force = false;

	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons &, TypeForces const &){
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){
	}

};

struct impulsive {

	static constexpr bool static_ions = false;
	static constexpr bool needs_force = false;	

	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons& ions, TypeForces const &){
		for(int i = 0; i != ions.num_atoms(); ++i)
			ions.coordinates()[i] += dt*ions.velocities()[i];
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){}

};


struct molecular_dynamics{

	static constexpr bool static_ions = false;
	static constexpr bool needs_force = true;	

	template <typename TypeIons, typename TypeForces>
	static auto acceleration(TypeIons& ions, TypeForces forces){

		for(int iatom = 0; iatom < ions.num_atoms(); iatom++) forces[iatom] /= ions.atoms()[iatom].mass();
		return forces;

	}
	
	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons& ions, TypeForces const & forces){
		solvers::velocity_verlet::propagate_positions(dt, acceleration(ions, forces), ions.velocities(), ions.coordinates());
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
