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

#ifndef PATH_HPP
#define PATH_HPP

#include <config.h>
#include <string>

namespace config {
  struct path {
    static std::string share(){ return SHARE_DIR + std::string("/") ; }
    static std::string unit_tests_data(){ return share() + std::string("unit_tests_data/"); }
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class config::path", "[path]") {
  SECTION("Share path"){
    REQUIRE(config::path::share() == SHARE_DIR + std::string("/"));
  }
}

#endif

#endif
