/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__GAUGE
#define INQ__PERTURBATIONS__GAUGE

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

enum class gauge {
  mixed,    // length gauge in non-periodic dimensions, velocity in periodic dimensions
  length,   // the electric field enters through the scalar potential or a phase in the orbitals, does not work for periodic dimensions
  velocity  // the electric field is applied through a uniform vector potential
};

}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_GAUGE_UNIT_TEST
#undef INQ_PERTURBATIONS_GAUGE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("perturbations::gauge", "[perturbations::gauge]") {

  auto gau = perturbations::gauge::velocity;
  
}

#endif
#endif
