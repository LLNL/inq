/* -*- indent-tabs-mode: t -*- */

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

#include <math/d3vector.hpp>
#include <utils/merge_optional.hpp>
#include <cassert>
#include <array>
#include <utils/merge_optional.hpp>

namespace input {

  class basis {

  public:

		static auto spacing(double arg_spacing){
			basis bs;
			bs.spacing_ = arg_spacing; 
			return bs;
		}

		static auto cutoff_energy(double arg_ecut){
			basis bs;
			bs.spacing_ = M_PI*sqrt(0.5/arg_ecut);
			return bs;
		}

		static auto spherical_grid(bool arg_sph_grid){
			basis bs;
			bs.spherical_grid_ = arg_sph_grid;
			return bs;
		}

		auto spacing() const {
			return spacing_.value();
		}

		auto spherical_grid() const {
			return spherical_grid_.value_or(false);
		}

		friend basis operator|(const basis & opt1, const basis & opt2){
			using utils::merge_optional;

			basis ropt;
			ropt.spacing_ = merge_optional(opt1.spacing_, opt2.spacing_);
			ropt.spherical_grid_ = merge_optional(opt1.spherical_grid_, opt2.spherical_grid_);
			return ropt;
		}

	private:
	
		basis(){
		}
		
		nonstd::optional<double> spacing_;
		nonstd::optional<bool> spherical_grid_;
		
	};
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class input::basis", "[basis]") {
  
  using namespace Catch::literals;

	SECTION("Spacing"){

		auto bi = input::basis::spacing(0.123);

		REQUIRE(bi.spacing() == 0.123_a);
		REQUIRE(not bi.spherical_grid());
				
	}
					
	
	SECTION("Cutoff energy"){

		auto bi = input::basis::cutoff_energy(493.48);

		REQUIRE(bi.spacing() == 0.1_a);
		
	}
			
	SECTION("Spherical grid"){

		auto bi = input::basis::cutoff_energy(493.48) | input::basis::spherical_grid(true);

		REQUIRE(bi.spacing() == 0.1_a);
		REQUIRE(bi.spherical_grid());
		
	}
			
}
#endif

    
#endif
