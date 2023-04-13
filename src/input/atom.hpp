/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__ATOM
#define INQ__INPUT__ATOM

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/vector3.hpp>

#include <pseudopod/element.hpp>
#include <input/species.hpp>
#include <vector>
#include <cmath>

namespace inq {
namespace input {

class atom {

public:
		
	atom(const input::species & arg_spec, const vector3<double> & arg_position):
		species_(arg_spec),
		position_(arg_position){
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

	friend auto operator==(atom const& self, atom const& other) {
		return
			    self.species_  == other.species_
			and self.position_ == other.position_
		;
	}
	friend auto operator!=(atom const& self, atom const& other) {return not(self == other);}

private:

	input::species species_;
	vector3<double> position_;

};
}
}

auto operator|(const inq::input::species & arg_spec, const inq::vector3<double> & arg_position){
	return inq::input::atom(arg_spec, arg_position);
}

auto operator|(const pseudo::element & arg_element, const inq::vector3<double> & arg_position){
	return inq::input::atom(arg_element, arg_position);
}

auto operator|(const std::string & arg_symbol, const inq::vector3<double> & arg_position){
	return inq::input::atom(pseudo::element(arg_symbol), arg_position);
}
#endif

#ifdef INQ_INPUT_ATOM_UNIT_TEST
#undef INQ_INPUT_ATOM_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace Catch::literals;

	SECTION("Constructor"){
		inq::input::atom at(pseudo::element("H"), inq::vector3<double>(1.0, 2.0, 3.0));

		CHECK(at.species().atomic_number() == 1);
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
	}
	
	SECTION("Species composition"){
		auto at = inq::input::species(pseudo::element("C")) | inq::vector3<double>(1.0, 2.0, 3.0);

		CHECK(at.species().symbol() == "C");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
	}
	
	SECTION("Species option composition"){
		
		auto at = pseudo::element("C") | inq::input::species::symbol("C1") | inq::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "C1");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
	}

	SECTION("Element composition"){
		
		auto at = pseudo::element("W") | inq::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "W");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
	}

	SECTION("String composition"){
		
		auto at = std::string("Xe") | inq::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "Xe");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
	}
	
	SECTION("Char * composition"){
		
		inq::input::atom at = "Tc" | inq::vector3<double>(1.0, 2.0, 3.0);
		
		CHECK(at.species().symbol() == "Tc");
		CHECK(at.position()[0] == 1.0_a);
		CHECK(at.position()[1] == 2.0_a);
		CHECK(at.position()[2] == 3.0_a);
	}

	SECTION("Equality") {

		inq::input::atom at1{pseudo::element("H"), inq::vector3<double>(1.0, 2.0, 3.0)};
		inq::input::atom at2{pseudo::element("H"), inq::vector3<double>(1.0, 2.0, 3.0)};

		CHECK( at1 == at2 );
	}
}
#endif
