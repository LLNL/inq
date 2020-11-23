/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__FORCES
#define INQ__HAMILTONIAN__FORCES

#include <systems/ions.hpp>
#include <systems/electrons.hpp>

/*
 Copyright (C) 2019 Xavier Andrade

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

namespace inq {
namespace hamiltonian {

auto calculate_forces(const systems::ions & ions, systems::electrons & electrons){
  math::array<math::vector3<double>, 1> forces(ions.geo().num_atoms(), {0.0, 0.0, 0.0});
  
  return forces;
}

}
}

#ifdef INQ_HAMILTONIAN_FORCES_UNIT_TEST
#undef INQ_HAMILTONIAN_FORCES_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::forces", "[forces]"){

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif

