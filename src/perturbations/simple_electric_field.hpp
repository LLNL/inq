/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__SIMPLE_ELECTRIC_FIELD
#define INQ__PERTURBATIONS__SIMPLE_ELECTRIC_FIELD

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

class simple_electric_field : public none {

	vector3<double> electric_field_;
	
public:

	simple_electric_field(vector3<double> const & electric_field):
		electric_field_(electric_field)
	{
	}
	
	auto has_uniform_electric_field() const {
		return true;
	}

	auto uniform_electric_field(double /*time*/) const {
		return electric_field_;
	}
	
	template<class OStream>
	friend OStream & operator<<(OStream & out, simple_electric_field const & self){
		return out;
	}
	
};
	
}
}
#endif

#ifdef INQ_PERTURBATIONS_SIMPLE_ELECTRIC_FIELD_UNIT_TEST
#undef INQ_PERTURBATIONS_SIMPLE_ELECTRIC_FIELD_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	perturbations::simple_electric_field efield{{1.0, 2.0, 3.0}};

	CHECK(efield.has_uniform_electric_field());
	CHECK(efield.uniform_electric_field(/*time = */ 0.0)    == vector3<double>{1.0, 2.0, 3.0});
	CHECK(efield.uniform_electric_field(/*time = */ 1000.0) == vector3<double>{1.0, 2.0, 3.0});
	
}
#endif
