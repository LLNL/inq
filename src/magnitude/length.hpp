/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__LENGTH
#define INQ__MAGNITUDE__LENGTH

/*
 Copyright (C) 2021 Xavier Andrade

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

#include <inq/quantity.hpp>

namespace inq {
namespace magnitude {

class length {
	
};

auto operator "" _b(long double val){
	return inq::quantity<length>::from_atomic_units(val);
}

auto operator "" _bohr(long double val){
	return inq::quantity<length>::from_atomic_units(val);
}

auto operator "" _angstrom(long double val){
	return inq::quantity<length>::from_atomic_units(1.88972612462938);
}

auto operator "" _Angstrom(long double val){
	return val*1.0_angstrom;
}

auto operator "" _A(long double val){
	return val*1.0_angstrom;
}

auto operator "" _nanometer(long double val){
	return val*10.0*1.0_angstrom;
}

auto operator "" _nm(long double val){
	return val*1.0_nanometer;
}

auto operator "" _picometer(long double val){
	return val/1000.0*1.0_nanometer;
}

auto operator "" _pm(long double val){
	return val*1.0_picometer;
}


}
}

#ifdef INQ_MAGNITUDE_LENGTH_UNIT_TEST
#undef INQ_MAGNITUDE_LENGTH_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("magnitude::quantity", "[magnitude::quantity]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace magnitude;

	auto l1 = 100.0_b;
	CHECK(l1.in_atomic_units() == 100.0);

  auto l2 = 100.0_A;
	CHECK(l2.in_atomic_units() == 188.9726124629_a);

  auto l3 = 23.333_pm;
  CHECK(l3.in_atomic_units() == 0.440929796659773_a);
	
}

#endif

#endif

