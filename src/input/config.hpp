/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__CONFIG
#define INQ__INPUT__CONFIG

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa

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

#include <cassert>
#include <magnitude/energy.hpp>

namespace inq {
namespace input {

struct config {

	int extra_states = 0;
	double excess_charge = 0.0;
	quantity<magnitude::energy> temperature = quantity<magnitude::energy>::zero();

};

}
}

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_CONFIG_UNIT_TEST
#undef INQ_INPUT_CONFIG_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("input::config", "[input::config]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
