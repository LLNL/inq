/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__BASIS
#define INQ__INPUT__BASIS

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

#include <math/vector3.hpp>
#include <utils/merge_optional.hpp>
#include <cassert>
#include <array>
#include <optional>
#include <utils/merge_optional.hpp>
#include <magnitude/energy.hpp>

namespace inq {
namespace input {

  class basis {

  public:

		static auto spacing(double arg_spacing){
			basis bs;
			bs.spacing_ = arg_spacing; 
			return bs;
		}

		static auto cutoff_energy(quantity<magnitude::energy> arg_ecut){
			basis bs;
			bs.spacing_ = M_PI*sqrt(0.5/arg_ecut.in_atomic_units());
			return bs;
		}

		static auto spherical_grid(bool arg_sph_grid){
			basis bs;
			bs.spherical_grid_ = arg_sph_grid;
			return bs;
		}

		static auto density_factor(double arg_factor){
			basis bs;
			bs.density_factor_ = arg_factor;
			return bs;
		}
		
		auto spacing_value() const {
			return spacing_.value();
		}

		auto spherical_grid_value() const {
			return spherical_grid_.value_or(false);
		}

		auto density_factor_value() const {
			return density_factor_.value_or(1.0);
		}

		static auto double_grid(){
			basis bs;
			bs.double_grid_ = true;
			return bs;
		}
		
		auto double_grid_value() const {
			return double_grid_.value_or(false);
		}
		
		friend basis operator|(const basis & opt1, const basis & opt2){
			using inq::utils::merge_optional;

			basis ropt;
			ropt.spacing_ = merge_optional(opt1.spacing_, opt2.spacing_);
			ropt.spherical_grid_ = merge_optional(opt1.spherical_grid_, opt2.spherical_grid_);
			ropt.density_factor_ = merge_optional(opt1.density_factor_, opt2.density_factor_);
			ropt.double_grid_	= merge_optional(opt1.double_grid_, opt2.double_grid_);		
			return ropt;
		}

		

	private:
	
		basis(){
		}
		
		std::optional<double> spacing_;
		std::optional<bool> spherical_grid_;
		std::optional<double> density_factor_;
		std::optional<bool> double_grid_;	

	};

}
}

#ifdef INQ_INPUT_BASIS_UNIT_TEST
#undef INQ_INPUT_BASIS_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class input::basis", "[basis]") {
  
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;

	SECTION("Spacing"){

		auto bi = input::basis::spacing(0.123);

		CHECK(bi.spacing_value() == 0.123_a);
		CHECK(not bi.spherical_grid_value());
				
	}
					
	
	SECTION("Cutoff energy"){

		auto bi = input::basis::cutoff_energy(493.48_Ha);

		CHECK(bi.spacing_value() == 0.1_a);
		
	}
			
	SECTION("Spherical grid"){

		auto bi = input::basis::cutoff_energy(493.48_Ha) | input::basis::spherical_grid(true);

		CHECK(bi.spacing_value() == 0.1_a);
		CHECK(bi.spherical_grid_value());
		
	}
			
}
#endif

    
#endif
