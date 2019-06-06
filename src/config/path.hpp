#ifndef PATH_HPP
#define PATH_HPP

namespace config {
  struct path {
    static std::string share(){ return SHARE_DIR; }
  };
}

#ifdef UNIT_TEST
#include <catch.hpp>

TEST_CASE("class config::path", "[path]") {
  SECTION("Share directory"){
    REQUIRE(config::path::share() == SHARE_DIR);
  }
}
 
#endif

#endif
