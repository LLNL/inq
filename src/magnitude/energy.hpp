/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__ENERGY
#define INQ__MAGNITUDE__ENERGY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/quantity.hpp>

namespace inq {
namespace magnitude {

struct energy;

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

static auto const Ha = inq::magnitude::operator""_Ha(1);

struct energy {

	static inq::quantity<energy> parse(std::string units){
		
		std::transform(units.begin(), units.end(), units.begin(), ::tolower);
		
		if(units == "hartree" or units == "hartrees" or units == "ha") {
			return 1.0_Ha;
		} else if (units == "electronvolt" or units == "electronvolts" or units == "ev"){
			return 1.0_eV;
		} else if (units == "rydberg" or units == "rydbergs" or units == "ry"){
			return 1.0_Ry;
		} else if (units == "kelvin" or units == "kelvins" or units == "k"){
			return 1.0_K;
		} else {
			throw std::runtime_error("inq error: unknown energy units '" + units + "'.");
		}
	}
	
};


}
}
#endif

#ifdef INQ_MAGNITUDE_ENERGY_UNIT_TEST
#undef INQ_MAGNITUDE_ENERGY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

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

	CHECK(energy::parse("hartree") == 1.0_Ha);	
	CHECK(energy::parse("Hartree") == 1.0_Ha);
	CHECK(energy::parse("HARTREE") == 1.0_Ha);
	CHECK(energy::parse("Ha") == 1.0_Ha);
	CHECK(energy::parse("electronvolt") == 1.0_eV);
	CHECK(energy::parse("eV") == 1.0_eV);
	CHECK(energy::parse("rydberg") == 1.0_Ry);	
	CHECK(energy::parse("Ry") == 1.0_Ry);
	CHECK(energy::parse("Kelvins") == 1.0_K);	
	CHECK(energy::parse("K") == 1.0_K);
	
	CHECK_THROWS(energy::parse("not_a_unit"));	
	
}
#endif

