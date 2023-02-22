/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__DYNAMICS
#define INQ__PERTURBATIONS__DYNAMICS

/*
 Copyright (C) 2022 Xavier Andrade

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
#include <basis/real_space.hpp>
#include <states/orbital_set.hpp>

namespace inq {
namespace perturbations {

enum class dynamics {
	none,
	polarization
};

}
}
#endif

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_DYNAMICS_UNIT_TEST
#undef INQ_PERTURBATIONS_DYNAMICS_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  [[maybe_unused]] auto dyn = perturbations::dynamics::none;

}

#endif
