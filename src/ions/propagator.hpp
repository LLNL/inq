/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__IONS__PROPAGATOR
#define INQ__IONS__PROPAGATOR

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
		for(int i = 0; i != ions.geo().num_atoms(); ++i)
			ions.geo().coordinates()[i] += dt*ions.geo().velocities()[i];
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){}

};


struct molecular_dynamics{

	static constexpr bool static_ions = false;
	static constexpr bool needs_force = true;	

	template <typename TypeIons, typename TypeForces>
	static auto acceleration(TypeIons& ions, TypeForces forces){

		for(int iatom = 0; iatom < ions.geo().num_atoms(); iatom++) forces[iatom] /= ions.geo().atoms()[iatom].mass();
		return forces;

	}
	
	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons& ions, TypeForces const & forces){
		solvers::velocity_verlet::propagate_positions(dt, acceleration(ions, forces), ions.geo().velocities(), ions.geo().coordinates());
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons & ions, TypeForces const & forces){
		solvers::velocity_verlet::propagate_velocities(dt, acceleration(ions, forces), ions.geo().velocities());
	}


};



}
}
}

#ifdef INQ_IONS_PROPAGATOR_UNIT_TEST
#undef INQ_IONS_PROPAGATOR_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("ions::propagator", "[ions::progagator]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
