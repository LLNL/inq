/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__NUM_STR
#define INQ__UTILS__NUM_STR

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cstdlib>

namespace inq {
namespace util {

auto num_to_str(long num){
	char numcstr[12]; 
	snprintf(numcstr, 11, "%010ld", num);
	return std::string(numcstr);				
}

template <typename Type>
auto str_to(std::string const & str);

template <>
auto str_to<double>(std::string const & str) {
  return std::atof(str.c_str());
}

template <>
auto str_to<int>(std::string const & str) {
  return std::atoi(str.c_str());
}

template <>
auto str_to<long>(std::string const & str) {
  return std::atol(str.c_str());
}

}
}
#endif

#ifdef INQ_UTILS_NUM_STR_UNIT_TEST
#undef INQ_UTILS_NUM_STR_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <math/complex.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

  SECTION("number to string"){
    
    CHECK(util::num_to_str(0)    == "0000000000");
    CHECK(util::num_to_str(1024) == "0000001024");
    CHECK(util::num_to_str(-333) == "-000000333");

  }

  SECTION("string to number"){

    CHECK(util::str_to<double>("0.3") == 0.3_a);
    CHECK(util::str_to<double>("-56750.25456") == -56750.25456_a);

    CHECK(util::str_to<int>("333") == 333);
    CHECK(util::str_to<int>("-45") == -45);

    CHECK(util::str_to<long>("46785454537460193") == 46785454537460193l);
    CHECK(util::str_to<long>("-467489690221") == -467489690221l);
  }
  

  
}

#endif

