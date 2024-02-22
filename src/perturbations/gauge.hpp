/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__GAUGE
#define INQ__PERTURBATIONS__GAUGE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <states/orbital_set.hpp>

namespace inq {
namespace perturbations {

enum class gauge {
  automatic,  // length gauge for non-periodic systems, velocity for periodic and partially-periodic systems
  length,     // the electric field enters through the scalar potential or a phase in the orbitals, does not work for periodic dimensions
  velocity    // the electric field is applied through a uniform vector potential
};

template<class OStream>
OStream & operator<<(OStream & out, gauge const & self){
	if(self == gauge::automatic)    out << "automatic";
	if(self == gauge::length)       out << "length";
	if(self == gauge::velocity)     out << "velocity";		
	return out;
}

template<class IStream>
IStream & operator>>(IStream & in, gauge & self){
	std::string readval;
	in >> readval;
	if(readval == "automatic"){
		self = gauge::automatic;
	} else if(readval == "length"){
		self = gauge::length;
	} else if(readval == "velocity"){
		self = gauge::velocity;
	} else {
		throw std::runtime_error("INQ error: Invalid gauge");
	}
	return in;
}

}
}
#endif

#ifdef INQ_PERTURBATIONS_GAUGE_UNIT_TEST
#undef INQ_PERTURBATIONS_GAUGE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace Catch::literals;
	using namespace inq;

	SECTION("Output"){
		{
			std::stringstream ss;
			std::string str;
			ss << perturbations::gauge::automatic;
			ss >> str;
			CHECK(str == "automatic");
		}

		{
			std::stringstream ss;
			std::string str;
			ss << perturbations::gauge::length;
			ss >> str;
			CHECK(str == "length");
		}

		{
			std::stringstream ss;
			std::string str;
			ss << perturbations::gauge::velocity;
			ss >> str;
			CHECK(str == "velocity");
		}
		
	}
	
	SECTION("Input"){
		{
			std::stringstream ss;
			perturbations::gauge sp;
			ss << "automatic";
			ss >> sp;
			CHECK(sp == perturbations::gauge::automatic);
		}

		{
			std::stringstream ss;
			perturbations::gauge sp;
			ss << "length";
			ss >> sp;
			CHECK(sp == perturbations::gauge::length);
		}
		
		{
			std::stringstream ss;
			perturbations::gauge sp;
			ss << "velocity";
			ss >> sp;
			CHECK(sp == perturbations::gauge::velocity);
		}
	}
	
	SECTION("Input/Output"){
		{
			std::stringstream ss;
			perturbations::gauge sp;
			ss << perturbations::gauge::automatic;
			ss >> sp;
			CHECK(sp == perturbations::gauge::automatic);
		}

		{
			std::stringstream ss;
			perturbations::gauge sp;
			ss << perturbations::gauge::length;
			ss >> sp;
			CHECK(sp == perturbations::gauge::length);
		}
		
		{
			std::stringstream ss;
			perturbations::gauge sp;
			ss << perturbations::gauge::velocity;
			ss >> sp;
			CHECK(sp == perturbations::gauge::velocity);
		}
	}
	
}
#endif
