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
	
}
#endif
