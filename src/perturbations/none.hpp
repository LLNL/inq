/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__NONE
#define INQ__PERTURBATIONS__NONE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>

namespace inq {
namespace perturbations {

class none {

public:

	template <typename DummyType>
	void zero_step(DummyType &) const {
	}
	
	auto has_uniform_electric_field() const {
		return false;
	}

	auto uniform_electric_field(double /*time*/) const {
		throw std::logic_error("This function should not be called");
		return vector3<double, cartesian>{0.0, 0.0, 0.0};
	}
	
	auto has_uniform_vector_potential() const {
		return false;
	}

	auto uniform_vector_potential(double /*time*/) const {
		throw std::logic_error("This function should not be called");
		return vector3<double, cartesian>{0.0, 0.0, 0.0};
	}
	
	template<typename PotentialType>
	void potential(const double time, PotentialType & potential) const {
	}

	template<class OStream>
	friend OStream & operator<<(OStream & out, none const & self){
		return out;
	}
	
};
	
}
}
#endif

#ifdef INQ_PERTURBATIONS_NONE_UNIT_TEST
#undef INQ_PERTURBATIONS_NONE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	perturbations::none nop;
	CHECK(not nop.has_uniform_electric_field());
}
#endif
