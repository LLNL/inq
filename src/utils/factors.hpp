/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__FACTORS
#define INQ__UTILS__FACTORS

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

#include <cassert>
#include <array>
#include <cmath>

#include <tinyformat/tinyformat.h>

namespace inq {
namespace utils {

  class factors {

  private:
   
    int number_;
    int current_;

    struct end_type {
    };
		
  public:
		
    factors(int number):
      number_(number),
      current_(1){
    }

    auto operator*() const {
      return current_;
    }

    void operator++(){
      while(current_ < number_) {
        current_++;
        if(number_%current_ == 0) break;
      }
    }

    static auto end(){
      return end_type{};
    }

    auto operator!=(end_type) const {
      return current_ != number_;
    }
		
  };

  class factors_reverse {

  private:
   
    int number_;
    int current_;

    struct end_type {
    };

  public:
		
    factors_reverse(int number):
      number_(number),
      current_(number){
    }

    auto operator*() const {
      return current_;
    }

    void operator++(){
      while(current_ > 1) {
        current_--;
        if(number_%current_ == 0) break;
      }
    }

    static auto end(){
      return end_type{};
    }

    auto operator!=(end_type) const {
      return current_ != 1;
    }
		
  };

}
}
#endif

#ifdef INQ_UTILS_FACTORS_UNIT_TEST
#undef INQ_UTILS_FACTORS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("class utils::factors", "[utils::factors]") {


  SECTION("Factors of 10"){
    inq::utils::factors fac(3);
    
    CHECK(*fac == 1);
    ++fac;
    CHECK(*fac == 3);
    CHECK(not (fac != fac.end()));
  }

  SECTION("Factors of 10"){
    inq::utils::factors fac(10);
    
    CHECK(*fac == 1);
    ++fac;
    CHECK(*fac == 2);
    ++fac;
    CHECK(fac != fac.end());  
    CHECK(*fac == 5);
    ++fac;
    CHECK(*fac == 10);
    CHECK(not (fac != fac.end()));
    
    ++fac; //check that it doesn't enter an infinite loop
    CHECK(*fac == 10);
  }

  SECTION("Factors of 12"){
    inq::utils::factors fac(12);
    
    CHECK(*fac == 1);
    ++fac;
    CHECK(*fac == 2);
    ++fac;
    CHECK(*fac == 3);
    ++fac;
    CHECK(*fac == 4);
    ++fac;
    CHECK(*fac == 6);
    ++fac;
    CHECK(*fac == 12);
    CHECK(not (fac != fac.end()));
    
    ++fac; //check that it doesn't enter an infinite loop
  }

	SECTION("Factors of 12 reversed"){
		inq::utils::factors_reverse fac(12);
    
    CHECK(*fac == 12);
    ++fac;
    CHECK(*fac == 6);
    ++fac;
    CHECK(*fac == 4);
    ++fac;
    CHECK(*fac == 3);
    ++fac;
    CHECK(*fac == 2);
    ++fac;
    CHECK(*fac == 1);
    CHECK(not (fac != fac.end()));
    
    ++fac; //check that it doesn't enter an infinite loop
  }


}
#endif
