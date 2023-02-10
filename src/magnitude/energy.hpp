/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__ENERGY
#define INQ__MAGNITUDE__ENERGY

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

class energy {
	
};

auto operator "" _ha(long double val){
	return inq::quantity<energy>::from_atomic_units(val);
}

auto operator "" _Ha(long double val){
	return inq::quantity<energy>::from_atomic_units(val);
}
	
auto operator "" _hartree(long double val){
	return inq::quantity<energy>::from_atomic_units(val);
}

auto operator "" _ev(long double val){
	return inq::quantity<energy>::from_atomic_units(0.0367493221756544*val);
}

auto operator "" _eV(long double val){
	return val*1.0_ev;
}
	
auto operator "" _electronvolt(long double val){
	return val*1.0_ev;
}

auto operator "" _ry(long double val){
	return inq::quantity<energy>::from_atomic_units(0.5*val);
}

auto operator "" _rydberg(long double val){
	return val*1.0_ry;
}

auto operator "" _Ry(long double val){
	return val*1.0_ry;
}
	
auto operator "" _K(long double val){
	return inq::quantity<energy>::from_atomic_units(3.16681156345556e-06*val);	
}
	
auto operator "" _kelvin(long double val){
	return val*1.0_K;
}
	
}
}
#endif

#ifdef INQ_MAGNITUDE_ENERGY_UNIT_TEST
#undef INQ_MAGNITUDE_ENERGY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("magnitude::energy", "[magnitude::energy]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace magnitude;

	{
		auto en = 100.0_ha;
		CHECK(en.in_atomic_units() == 100.0);
	}

	{
		auto en = 100.0_Ha;
		CHECK(en.in_atomic_units() == 100.0);
	}

	{
		auto en = 100.0_hartree;
		CHECK(en.in_atomic_units() == 100.0);
	}

	{
		auto en = 33.3_rydberg;
		CHECK(en.in_atomic_units() == 16.65);
	}

 	{
		auto en = 33.3_ry;
		CHECK(en.in_atomic_units() == 16.65);
	}

	{
		auto en = 33.3_Ry;
		CHECK(en.in_atomic_units() == 16.65);
	}
		
	{
		auto en = 300.0_K;
		CHECK(en.in_atomic_units() == 0.0009500435_a);
	}
	
	{
		auto en = 300.0_kelvin;
		CHECK(en.in_atomic_units() == 0.0009500435_a);
	}

	{
		auto en = 0.5_hartree + 300.0_kelvin;
		CHECK(en.in_atomic_units() == 0.5009500435_a);
	}
	
}
#endif

