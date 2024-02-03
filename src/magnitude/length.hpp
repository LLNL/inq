/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__LENGTH
#define INQ__MAGNITUDE__LENGTH

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/quantity.hpp>
#include <utils/lowercase.hpp>

namespace inq {
namespace magnitude {

struct length;

auto operator "" _b(long double val){
	return inq::quantity<length>::from_atomic_units(val);
}

auto operator "" _bohr(long double val){
	return inq::quantity<length>::from_atomic_units(val);
}

auto operator "" _angstrom(long double val){
	return inq::quantity<length>::from_atomic_units(1.88972612462938*val);
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

static auto const bohr     = inq::magnitude::operator""_bohr(1);

static auto const angstrom = inq::magnitude::operator""_angstrom(1); // Å, AA, Angstrom
static auto const A        = inq::magnitude::operator""_A(1);        // Å, AA, Angstrom

struct length {
	static inq::quantity<length> parse(std::string units){

		units = utils::lowercase(units);
		
		if(units == "bohr" or units == "bohrs" or units == "b") {
			return 1.0_b;
		} else if (units == "angstrom" or units == "angstroms" or units == "a"){
			return 1.0_A;
		} else if (units == "nanometer" or units == "nanometers" or units == "nm"){
			return 1.0_nm;
		} else if (units == "picometer" or units == "picometers" or units == "pm"){
			return 1.0_pm;
		} else {
			throw std::runtime_error("inq error: unknown length units '" + units + "'.");
		}
	}
	
};

}

}

#endif

#ifdef INQ_MAGNITUDE_LENGTH_UNIT_TEST
#undef INQ_MAGNITUDE_LENGTH_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using namespace magnitude;

  {
    auto le = 100.0_b;
    CHECK(le.in_atomic_units() == 100.0);
  }

  {
    auto le = 100.0_bohr;
    CHECK(le.in_atomic_units() == 100.0);
  }
  
  {
    auto le = 100.0_A;
    CHECK(le.in_atomic_units() == 188.9726124629_a);
  }
  
  {
    auto le = 100.0_angstrom;
    CHECK(le.in_atomic_units() == 188.9726124629_a);
  }

  {
    auto le = 874.21_nanometer;
    CHECK(le.in_atomic_units() == 16520.1747541225_a);
  }

  {
    auto le = 874.21_nm;
    CHECK(le.in_atomic_units() == 16520.1747541225_a);
  }
   
  {
    auto le = 23.333_picometer;
    CHECK(le.in_atomic_units() == 0.440929796659773_a);
  }
  
  {
    auto le = 23.333_pm;
    CHECK(le.in_atomic_units() == 0.440929796659773_a);
  }

	CHECK(length::parse("bohr") == 1.0_b);	
	CHECK(length::parse("Bohr") == 1.0_b);
	CHECK(length::parse("BOHR") == 1.0_b);
	CHECK(length::parse("a") == 1.0_A);
	CHECK(length::parse("angstrom") == 1.0_A);
	CHECK(length::parse("angSTROMS") == 1.0_A);		
	CHECK(length::parse("PICOMETER") == 1.0_pm);	
	CHECK(length::parse("picometers") == 1.0_pm);
	CHECK(length::parse("pm") == 1.0_pm);	
	CHECK(length::parse("NANOmeter") == 1.0_nm);	
	
	CHECK_THROWS(length::parse("not_a_unit") == 1.0_b);	
	
}
#endif
