/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__STATES__SPIN_CONFIG
#define INQ__STATES__SPIN_CONFIG

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/complex.hpp>
#include <gpu/array.hpp>
#include <operations/sum.hpp>
#include <parallel/partition.hpp>

namespace inq {
namespace states {

enum class spin_config { UNPOLARIZED, POLARIZED, NON_COLLINEAR };

template<class OStream>
OStream & operator<<(OStream & out, spin_config const & self){
	
	if(self == states::spin_config::UNPOLARIZED){
		out << "unpolarized";
	} else if(self == states::spin_config::POLARIZED){
		out << "polarized";
	} else if(self == states::spin_config::NON_COLLINEAR){
		out << "non_collinear";
	}
	
	return out;
}

template<class IStream>
IStream & operator>>(IStream & in, spin_config & self){

	std::string readval;
	in >> readval;

	if(readval == "unpolarized"){
		self = states::spin_config::UNPOLARIZED;
	} else if(readval == "polarized"){
		self = states::spin_config::POLARIZED;
	} else if(readval == "non_collinear"){
		self = states::spin_config::NON_COLLINEAR;
	} else {
		throw std::runtime_error("INQ error: Invalid spin configuration");
	}
	
	return in;
}

}
}
#endif

#ifdef INQ_STATES_SPIN_CONFIG_UNIT_TEST
#undef INQ_STATES_SPIN_CONFIG_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace Catch::literals;
	using namespace inq;


	SECTION("Output"){
		{
			std::stringstream ss;
			std::string str;
			ss << states::spin_config::UNPOLARIZED;
			ss >> str;
			CHECK(str == "unpolarized");
		}

		{
			std::stringstream ss;
			std::string str;
			ss << states::spin_config::POLARIZED;
			ss >> str;
			CHECK(str == "polarized");
		}

		{
			std::stringstream ss;
			std::string str;
			ss << states::spin_config::NON_COLLINEAR;
			ss >> str;
			CHECK(str == "non_collinear");
		}
		
	}
	
	SECTION("Input"){
		{
			std::stringstream ss;
			states::spin_config sp;
			ss << "unpolarized";
			ss >> sp;
			CHECK(sp == states::spin_config::UNPOLARIZED);
		}

		{
			std::stringstream ss;
			states::spin_config sp;
			ss << "polarized";
			ss >> sp;
			CHECK(sp == states::spin_config::POLARIZED);
		}
		
		{
			std::stringstream ss;
			states::spin_config sp;
			ss << "non_collinear";
			ss >> sp;
			CHECK(sp == states::spin_config::NON_COLLINEAR);
		}
	}
	
	SECTION("Input/Output"){
		{
			std::stringstream ss;
			states::spin_config sp;
			ss << states::spin_config::UNPOLARIZED;
			ss >> sp;
			CHECK(sp == states::spin_config::UNPOLARIZED);
		}

		{
			std::stringstream ss;
			states::spin_config sp;
			ss << states::spin_config::POLARIZED;
			ss >> sp;
			CHECK(sp == states::spin_config::POLARIZED);
		}
		
		{
			std::stringstream ss;
			states::spin_config sp;
			ss << states::spin_config::NON_COLLINEAR;
			ss >> sp;
			CHECK(sp == states::spin_config::NON_COLLINEAR);
		}
	}
	
}
#endif
