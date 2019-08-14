/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__CELL
#define INPUT__CELL

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

#include <math/d3vector.hpp>
#include <cassert>
#include <array>

namespace input {

  class cell {

  public:

		static auto cubic(double aa){
			return cell(math::d3vector(aa, 0.0, 0.0), math::d3vector(0.0, aa, 0.0), math::d3vector(0.0, 0.0, aa));
		}

		static auto cubic(double aa, double bb, double cc){
			return cell(math::d3vector(aa, 0.0, 0.0), math::d3vector(0.0, bb, 0.0), math::d3vector(0.0, 0.0, cc));
		}

		auto & operator[](const int ii) const {
			return lattice_vectors_[ii].value();
		}
		
	private:

		cell(const math::d3vector & a0, const math::d3vector & a1, const math::d3vector & a2){
			lattice_vectors_[0] = a0;
			lattice_vectors_[1] = a1;
			lattice_vectors_[2] = a2;
		}

		std::array<std::optional<math::d3vector>, 3> lattice_vectors_;
		
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class input::cell", "[cell]") {
  
  using namespace Catch::literals;

	SECTION("Cubic"){

		auto ci = input::cell::cubic(10.2);

		REQUIRE(ci[0][0] == 10.2_a);
		REQUIRE(ci[0][1] == 0.0_a);
		REQUIRE(ci[0][2] == 0.0_a);
		REQUIRE(ci[1][0] == 0.0_a);
		REQUIRE(ci[1][1] == 10.2_a);
		REQUIRE(ci[1][2] == 0.0_a);
		REQUIRE(ci[2][0] == 0.0_a);
		REQUIRE(ci[2][1] == 0.0_a);
		REQUIRE(ci[2][2] == 10.2_a);
		
	}
	
	SECTION("Parallelepipedic"){

		auto ci = input::cell::cubic(10.2, 5.7, 8.3);

		REQUIRE(ci[0][0] == 10.2_a);
		REQUIRE(ci[0][1] == 0.0_a);
		REQUIRE(ci[0][2] == 0.0_a);
		REQUIRE(ci[1][0] == 0.0_a);
		REQUIRE(ci[1][1] == 5.7_a);
		REQUIRE(ci[1][2] == 0.0_a);
		REQUIRE(ci[2][0] == 0.0_a);
		REQUIRE(ci[2][1] == 0.0_a);
		REQUIRE(ci[2][2] == 8.3_a);
		
	}
  
}
#endif

    
#endif
