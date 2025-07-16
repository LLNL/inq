/* -*- indent-tabs-mode: t -*- */
#ifndef INQ__CONFIG__PATH
#define INQ__CONFIG__PATH

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <string>

namespace inq {
namespace config {
namespace path {

std::string share(){
	auto env = std::getenv("INQ_SHARE_PATH");
	if(env == NULL) return SHARE_DIR + std::string("/");
	return env + std::string("/");
}

std::string unit_tests_data(){
	return share() + std::string("unit_tests_data/");
}

}
}
}
#endif

#ifdef INQ_CONFIG_PATH_UNIT_TEST
#undef INQ_CONFIG_PATH_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	//Do not run if the variable is set
	if(std::getenv("INQ_SHARE_PATH") != NULL) return;
	
  SECTION("Share path"){
    CHECK(inq::config::path::share() == SHARE_DIR + std::string("/"));
  }

	SECTION("Unit test path"){
    CHECK(inq::config::path::unit_tests_data() == SHARE_DIR + std::string("/unit_tests_data/"));
  }
	
	SECTION("Path from environment variable"){
		REQUIRE(setenv("INQ_SHARE_PATH", "/basura", 1) == 0);
    CHECK(inq::config::path::share() == std::string("/basura/"));
		REQUIRE(unsetenv("INQ_SHARE_PATH") == 0);
  }

}

#endif
