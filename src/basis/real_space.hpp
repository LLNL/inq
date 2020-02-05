/* -*- indent-tabs-mode: t -*- */

#ifndef BASIS_REAL_SPACE
#define BASIS_REAL_SPACE

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

#include <math/vec3d.hpp>
#include <ions/unitcell.hpp>
#include "grid.hpp"
#include <cassert>
#include <array>
#include <input/basis.hpp>
#include <gpu/run.hpp>

namespace basis {

  class real_space : public grid {

  public:
		
    real_space(const ions::UnitCell & cell, const input::basis & basis_input, comm_type comm = MPI_COMM_SELF):
			grid(cell, calculate_dimensions(cell, basis_input), basis_input.spherical_grid(), cell.periodic_dimensions(), comm){
    }

		real_space(const grid & grid_basis):
			grid(grid_basis){
    }

		GPU_FUNCTION math::vec3d rvector(const int ix, const int iy, const int iz) const {
			math::vec3d rr{ix*rspacing()[0], iy*rspacing()[1], iz*rspacing()[2]};
			for(int idir = 0; idir < 3; idir++) {
				if(rr[idir] >= 0.5*rlength()[idir]) rr[idir] -= rlength()[idir];
			}
			return rr;
		}
		
		template <class int_array>
		GPU_FUNCTION math::vec3d rvector(const int_array & indices) const {
			return rvector(indices[0], indices[1], indices[2]);
		}
		
		double r2(const int ix, const int iy, const int iz) const {
			return norm(rvector(ix, iy, iz));
		}

		friend auto operator==(const real_space & rs1, const real_space & rs2){
			bool equal = rs1.nr_[0] == rs2.nr_[0] and rs1.nr_[1] == rs2.nr_[1] and rs1.nr_[2] == rs2.nr_[2];
			equal = equal and rs1.rspacing()[0] == rs2.rspacing()[0];
			equal = equal and rs1.rspacing()[1] == rs2.rspacing()[1];
			equal = equal and rs1.rspacing()[2] == rs2.rspacing()[2];
			return equal;
		}

		auto enlarge(int factor) const {
			return real_space(grid(cell_.enlarge(factor), {factor*nr_[0], factor*nr_[1], factor*nr_[2]}, spherical_g_grid_, periodic_dimensions_, this->dist().comm()));
		}

		auto volume_element() const {
			return rspacing_[0]*rspacing_[1]*rspacing_[2];
		}
		
	private:

		static std::array<int, 3> calculate_dimensions(const ions::UnitCell & cell, const input::basis & basis_input){
			std::array<int, 3> nr;
			double spacing = basis_input.spacing();
			
			// make the spacing conmensurate with the grid
			// OPTIMIZATION: we can select a good size here for the FFT
			for(int idir = 0; idir < 3; idir++){
				double rlength = length(cell[idir]);
				nr[idir] = round(rlength/spacing);
			}
			
			return nr;
		}
		
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::real_space", "[real_space]") {
  
  using namespace Catch::literals;
  using math::vec3d;
  
  {
    
    SECTION("Cubic cell"){

      ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 10.0, 0.0), vec3d(0.0, 0.0, 10.0));

      double ecut = 20.0;
      
      basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

      REQUIRE(rs.size() == 8000);
      
      REQUIRE(rs.rspacing()[0] == 0.5_a);
      REQUIRE(rs.rspacing()[1] == 0.5_a);
      REQUIRE(rs.rspacing()[2] == 0.5_a);
      
      REQUIRE(rs.sizes()[0] == 20);
      REQUIRE(rs.sizes()[1] == 20);
      REQUIRE(rs.sizes()[2] == 20);

    }

    SECTION("Parallelepipedic cell"){

      ions::UnitCell cell(vec3d(77.7, 0.0, 0.0), vec3d(0.0, 14.14, 0.0), vec3d(0.0, 0.0, 23.25));

      double ecut = 37.9423091;
      
      basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

      REQUIRE(rs.size() == 536640);
	    
      REQUIRE(rs.rspacing()[0] == 0.3613953488_a);
      REQUIRE(rs.rspacing()[1] == 0.3625641026_a);
      REQUIRE(rs.rspacing()[2] == 0.36328125_a);
      
      REQUIRE(rs.sizes()[0] == 215);
      REQUIRE(rs.sizes()[1] == 39);
			REQUIRE(rs.sizes()[2] == 64);

			auto rs3x = rs.enlarge(3);
			
      REQUIRE(rs3x.rspacing()[0] == 0.3613953488_a);
      REQUIRE(rs3x.rspacing()[1] == 0.3625641026_a);
      REQUIRE(rs3x.rspacing()[2] == 0.36328125_a);
      
      REQUIRE(rs3x.sizes()[0] == 3*215);
      REQUIRE(rs3x.sizes()[1] == 3*39);
			REQUIRE(rs3x.sizes()[2] == 3*64);
			
    }

  }
}
#endif

    
#endif
