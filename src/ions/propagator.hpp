/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__IONS__PROPAGATOR
#define INQ__IONS__PROPAGATOR

namespace inq {
namespace ions {
namespace propagator {

struct fixed {

	static constexpr bool static_ions = true;

	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons &, TypeForces const &){
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){
	}

};

struct impulsive {

	static constexpr bool static_ions = false;

	template <typename TypeIons, typename TypeForces>
	static void propagate_positions(double dt, TypeIons& ions, TypeForces const &){
		for(int i = 0; i != ions.geo().num_atoms(); ++i)
			ions.geo().coordinates()[i] += dt*ions.velocities()[i];
	}

	template <typename TypeIons, typename TypeForces>
	static void propagate_velocities(double dt, TypeIons &, TypeForces const &){}

};

}
}
}

#ifdef INQ_IONS_PROPAGATOR_UNIT_TEST
#undef INQ_IONS_PROPAGATOR_UNIT_TEST

#endif

#endif

