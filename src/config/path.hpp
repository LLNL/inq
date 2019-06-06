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
#include <catch.hpp>

TEST_CASE("class config::path", "[path]") {
  SECTION("Share path"){
    REQUIRE(config::path::share() == SHARE_DIR + std::string("/"));
  }
}

#endif

#endif
