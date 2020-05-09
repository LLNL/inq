/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__RT
#define INPUT__RT

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <utils/merge_optional.hpp>

#include <nonstd/optional.hpp>
#include <cassert>

namespace input {

  class rt {

  public:

    rt(){
		}

    static auto dt(double dt) {
      rt solver;
      solver.dt_ = dt;
      return solver;
    }

    auto dt() const {
      return dt_.value_or(0.01);
    }

		auto static num_steps(double etol) {
			rt solver;
      solver.num_steps_ = etol;
      return solver;
    }
				
 		auto num_steps() const {
			return num_steps_.value_or(100);
		}
        
    friend auto operator|(const rt & solver1, const rt & solver2){
			using utils::merge_optional;

			rt rsolver;
			rsolver.dt_	= merge_optional(solver1.dt_, solver2.dt_);
			rsolver.num_steps_	= merge_optional(solver1.num_steps_, solver2.num_steps_);
			return rsolver;
		}
    
  private:

    nonstd::optional<double> dt_;
		nonstd::optional<int> num_steps_;
		
  };
    
}

////////////////////////////////////////////////////////

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class input::rt", "[input::rt]") {
  
  using namespace Catch::literals;

	SECTION("Defaults"){

    input::rt solver;

    CHECK(solver.dt() == 0.01_a);
    CHECK(solver.num_steps() == 100);		
    
  }

  SECTION("Composition"){

    auto solver = input::rt::num_steps(1000) | input::rt::dt(0.05);
    
    CHECK(solver.num_steps() == 1000);
    CHECK(solver.dt() == 0.05_a);
    
  }

}

#endif
   
#endif
