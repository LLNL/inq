/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__ATOM
#define INQ__INPUT__ATOM

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

#include <pseudopod/element.hpp>
#include <input/species.hpp>
#include <vector>
#include <cmath>

namespace inq {
namespace input {

class atom {

public:
		
	atom(const input::species & arg_spec, const math::vec3d & arg_position):
		species_(arg_spec),
		position_(arg_position){
	}

	const auto & species() const {
		return species_;
	}

	const auto & position() const {
		return position_;
	}
		
private:

	input::species species_;
	math::vec3d position_;
		
};
}

auto operator|(const input::species & arg_spec, const math::vec3d & arg_position){
	return input::atom(arg_spec, arg_position);
}

auto operator|(const pseudo::element & arg_element, const math::vec3d & arg_position){
	return input::atom(arg_element, arg_position);
}

auto operator|(const std::string & arg_symbol, const math::vec3d & arg_position){
	return input::atom(pseudo::element(arg_symbol), arg_position);
}

}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class ions::atom", "[input::atom]") {
  
	using namespace inq;
	using namespace Catch::literals;

	SECTION("Constructor"){
		input::atom at(pseudo::element("H"), math::vec3d(1.0, 2.0, 3.0));

		CHECK(at.species().atomic_number() == 1);
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
		
	}
	
	SECTION("Species composition"){
		input::atom at = input::species(pseudo::element("C")) | math::vec3d(1.0, 2.0, 3.0);

		CHECK(at.species().symbol() == "C");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);

	}
	
	SECTION("Species option composition"){
		
		input::atom at = pseudo::element("C") | input::species::symbol("C1") | math::vec3d(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "C1");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);

	}

	SECTION("Element composition"){
		
		input::atom at = pseudo::element("W") | math::vec3d(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "W");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);

	}

	SECTION("String composition"){
		
		input::atom at = std::string("Xe") | math::vec3d(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "Xe");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
		
	}
	
	SECTION("Char * composition"){
		
		input::atom at = "Tc" | math::vec3d(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "Tc");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
		
	}
	
}


#endif

#endif
