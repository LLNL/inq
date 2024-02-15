/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__PROPAGATOR
#define INQ__IONS__PROPAGATOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <options/real_time.hpp>
#include <solvers/velocity_verlet.hpp>
#include <variant>

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
	void propagate_positions(double dt, TypeIons &, TypeForces const &) const {
	}

	template <typename TypeIons, typename TypeForces>
	void propagate_velocities(double dt, TypeIons &, TypeForces const &) const {
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
	void propagate_positions(double dt, TypeIons& ions, TypeForces const &) const {
		for(int i = 0; i != ions.size(); ++i)
			ions.positions()[i] += dt*ions.velocities()[i];
	}

	template <typename TypeIons, typename TypeForces>
	void propagate_velocities(double dt, TypeIons &, TypeForces const &) const {}

};


struct molecular_dynamics{

	constexpr bool static_ions() const {
		return false;
	}

	constexpr bool needs_force() const {
		return true;
	}

	template <typename TypeIons, typename TypeForces>
	auto acceleration(TypeIons& ions, TypeForces forces) const {

		for(int iatom = 0; iatom < ions.size(); iatom++) forces[iatom] /= ions.atoms()[iatom].mass();
		return forces;

	}
	
	template <typename TypeIons, typename TypeForces>
	void propagate_positions(double dt, TypeIons& ions, TypeForces const & forces) const {
		solvers::velocity_verlet::propagate_positions(dt, acceleration(ions, forces), ions.velocities(), ions.positions());
	}

	template <typename TypeIons, typename TypeForces>
	void propagate_velocities(double dt, TypeIons & ions, TypeForces const & forces) const {
		solvers::velocity_verlet::propagate_velocities(dt, acceleration(ions, forces), ions.velocities());
	}
};


class runtime {

	std::variant<fixed, impulsive, molecular_dynamics> var_;

public:
	
	runtime(options::real_time::ion_dynamics arg_dynamics) {
		switch (arg_dynamics) {
    case options::real_time::ion_dynamics::STATIC:
			var_ = ions::propagator::fixed{};
			break;
    case options::real_time::ion_dynamics::IMPULSIVE:
			var_ = ions::propagator::impulsive{};
			break;
    case options::real_time::ion_dynamics::EHRENFEST:
			var_ = ions::propagator::molecular_dynamics{};
			break;
		}
	}

	constexpr bool static_ions() const {
		return std::visit([&](auto ip) { return ip.static_ions(); }, var_);
	}
	
	constexpr bool needs_force() const {
		return std::visit([&](auto ip) { return ip.needs_force(); }, var_);
	}

	template <typename TypeIons, typename TypeForces>
	void propagate_positions(double dt, TypeIons& ions, TypeForces const & forces) const {
		return std::visit([&](auto ip) { return ip.propagate_positions(dt, ions, forces); }, var_);
	}
	
	template <typename TypeIons, typename TypeForces>
	void propagate_velocities(double dt, TypeIons& ions, TypeForces const & forces) const {
		return std::visit([&](auto ip) { return ip.propagate_velocities(dt, ions, forces); }, var_);
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
