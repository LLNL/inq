/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__COORD
#define INQ__INPUT__COORD

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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

#include <math/vector3.hpp>
#include <inq/quantity.hpp>
#include <magnitude/length.hpp>

namespace inq {
namespace input {

using coord = math::vector3<autocast_quantity<magnitude::length>>;

}
}

#ifdef INQ_INPUT_COORD_UNIT_TEST
#undef INQ_INPUT_COORD_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("input::coord", "[input::coord]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
