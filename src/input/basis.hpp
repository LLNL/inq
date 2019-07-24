/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef INPUT__BASIS
#define INPUT__BASIS

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

#include "../math/d3vector.hpp"
#include <cassert>
#include <array>

namespace input {

  class basis {

  public:

		static auto spacing(double arg_spacing){
			return basis(arg_spacing);
		}

		static auto cutoff_energy(double arg_ecut){
			return basis(M_PI*sqrt(0.5/arg_ecut));
		}

		double get_spacing() const {
			return spacing_;
		}
		
	private:

		basis(double arg_spacing):
			spacing_(arg_spacing){
		}
		
		double spacing_;
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class input::basis", "[basis]") {
  
  using namespace Catch::literals;

	SECTION("Spacing"){

		auto bi = input::basis::spacing(0.123);

		REQUIRE(bi.get_spacing() == 0.123_a);
	}
					
	
	SECTION("Cutoff energy"){

		auto bi = input::basis::cutoff_energy(493.48);

		REQUIRE(bi.get_spacing() == 0.1_a);
	}
			
  
}
#endif

    
#endif
