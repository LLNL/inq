/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__KPOINTS
#define INQ__INPUT__KPOINTS

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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

#include <math/vector3.hpp>

#include <vector>
#include <cmath>

namespace inq {
namespace input {

class kpoints {

public:

  static auto gamma(){
    return kpoints({1, 1, 1}, {0.0, 0.0, 0.0});
  }
  
  static auto grid(math::vector3<int> const & dims, math::vector3<double> const & shifts = {0.0, 0.0, 0.0}){
    return kpoints(dims, shifts);
  }

  auto & dims() const {
    return dims_;
  }

  auto & shifts() const {
    return shifts_;
  }   

  auto num() const {
    return product(dims_);
  }
  
private:

  kpoints(math::vector3<int> const & dims, math::vector3<double> const & shifts):
    dims_(dims),
    shifts_(shifts)
  {
  }
  
  math::vector3<int> dims_;
  math::vector3<double> shifts_;

};
}
}

#ifdef INQ_INPUT_KPOINTS_UNIT_TEST
#undef INQ_INPUT_KPOINTS_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("class ions::kpoints", "[inq::input::kpoints]") {

  using namespace inq;
	using namespace Catch::literals;

	SECTION("Gamma - no arguments"){
		auto kpts = input::kpoints::gamma();

    CHECK(kpts.dims()[0] == 1);
    CHECK(kpts.dims()[1] == 1);
    CHECK(kpts.dims()[2] == 1);
    CHECK(kpts.shifts()[0] == 0.0_a);
    CHECK(kpts.shifts()[1] == 0.0_a);
    CHECK(kpts.shifts()[2] == 0.0_a);
    
	}
  
	SECTION("Grid - one argument"){
		auto kpts = input::kpoints::grid({10, 9, 8});

    CHECK(kpts.dims()[0] == 10);
    CHECK(kpts.dims()[1] == 9);
    CHECK(kpts.dims()[2] == 8);
    CHECK(kpts.shifts()[0] == 0.0_a);
    CHECK(kpts.shifts()[1] == 0.0_a);
    CHECK(kpts.shifts()[2] == 0.0_a);
    
	}
	
	SECTION("Grid - two arguments"){
		auto kpts = input::kpoints::grid({10, 9, 8}, {0.2, 0.1, 0.5});

    CHECK(kpts.dims()[0] == 10);
    CHECK(kpts.dims()[1] == 9);
    CHECK(kpts.dims()[2] == 8);
    CHECK(kpts.shifts()[0] == 0.2_a);
    CHECK(kpts.shifts()[1] == 0.1_a);
    CHECK(kpts.shifts()[2] == 0.5_a);
    
	}
}


#endif

#endif
