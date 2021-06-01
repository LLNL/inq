/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__ATOM
#define INQ__INPUT__ATOM

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

#include <pseudopod/element.hpp>
#include <input/coord.hpp>
#include <input/species.hpp>
#include <vector>
#include <cmath>

namespace inq {
namespace input {

class atom {

public:
		
	atom(const input::species & arg_spec, const math::vector3<double> & arg_position):
		species_(arg_spec),
		position_(arg_position){
	}

	atom(const input::species & arg_spec, const inq::input::coord & arg_coord):
		species_(arg_spec),
		position_({arg_coord[0], arg_coord[1], arg_coord[2]}){
	}
	
	const auto & species() const {
		return species_;
	}

	const auto & position() const {
		return position_;
	}

	auto & position(){
		return position_;
	}

private:

	input::species species_;
	math::vector3<double> position_;

};
}
}

auto operator|(const inq::input::species & arg_spec, const inq::math::vector3<double> & arg_position){
	return inq::input::atom(arg_spec, arg_position);
}

auto operator|(const pseudo::element & arg_element, const inq::math::vector3<double> & arg_position){
	return inq::input::atom(arg_element, arg_position);
}

auto operator|(const std::string & arg_symbol, const inq::math::vector3<double> & arg_position){
	return inq::input::atom(pseudo::element(arg_symbol), arg_position);
}

auto operator|(const inq::input::species & arg_spec, const inq::input::coord & arg_position){
	return inq::input::atom(arg_spec, arg_position);
}

auto operator|(const pseudo::element & arg_element, const inq::input::coord & arg_position){
	return inq::input::atom(arg_element, arg_position);
}

auto operator|(const std::string & arg_symbol, const inq::input::coord & arg_position){
	return inq::input::atom(pseudo::element(arg_symbol), arg_position);
}

#ifdef INQ_INPUT_ATOM_UNIT_TEST
#undef INQ_INPUT_ATOM_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("class ions::atom", "[inq::input::atom]") {
  
	using namespace Catch::literals;

	SECTION("Constructor"){
		inq::input::atom at(pseudo::element("H"), inq::math::vector3<double>(1.0, 2.0, 3.0));

		CHECK(at.species().atomic_number() == 1);
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
		
	}
	
	SECTION("Species composition"){
		auto at = inq::input::species(pseudo::element("C")) | inq::math::vector3<double>(1.0, 2.0, 3.0);

		CHECK(at.species().symbol() == "C");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);

	}
	
	SECTION("Species option composition"){
		
		auto at = pseudo::element("C") | inq::input::species::symbol("C1") | inq::math::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "C1");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);

	}

	SECTION("Element composition"){
		
		auto at = pseudo::element("W") | inq::math::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "W");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);

	}

	SECTION("String composition"){
		
		auto at = std::string("Xe") | inq::math::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "Xe");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
		
	}
	
	SECTION("Char * composition"){
		
		inq::input::atom at = "Tc" | inq::math::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "Tc");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
		
	}
	
}


#endif

#endif
