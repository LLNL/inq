/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__NONE
#define INQ__PERTURBATIONS__NONE

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>

namespace inq {
namespace perturbations {

class none {

public:

	auto has_uniform_electric_field() const {
		return false;
	}

	auto uniform_electric_field(double time) const {
		return math::vector3<double, math::cartesian>{0.0, 0.0, 0.0};
	}

};
	
}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_NONE_UNIT_TEST
#undef INQ_PERTURBATIONS_NONE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE("perturbations::none", "[perturbations::none]") {
	perturbations::none nop;
	CHECK(not nop.has_uniform_electric_field());
}

#endif
#endif
