#ifndef PLANE_WAVE_HPP
#define PLANE_WAVE_HPP

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

namespace basis {

  class plane_wave {

  public:
		
		plane_wave(ions::UnitCell& cell, std::array<int, 3> nr) : nr_{nr}{
			for(int idir = 0; idir < 3; idir++){
				rlength_[idir] = length(cell[idir]);
				ng_[idir] = nr_[idir];
				rspacing_[idir] = rlength_[idir]/nr_[idir];
				glength_[idir] = M_PI/rspacing_[idir];
				gspacing_[idir] = glength_[idir]/ng_[idir];
			}
		}

    plane_wave(ions::UnitCell & cell, const double & ecut) {
      ecut_ = ecut;
      rspacing_ = math::d3vector(M_PI*sqrt(0.5/ecut));
			
      // make the spacing conmensurate with the grid
      // OPTIMIZATION: we can select a good size here for the FFT
      for(int idir = 0; idir < 3; idir++){
				rlength_[idir] = length(cell[idir]);
				
				nr_[idir] = round(rlength_[idir]/rspacing_[idir]);
				ng_[idir] = nr_[idir];
				
				rspacing_[idir] = rlength_[idir]/nr_[idir];
				glength_[idir] = M_PI/rspacing_[idir];
				gspacing_[idir] = glength_[idir]/ng_[idir];
      }
    }

    const double & ecut() const {
      return ecut_;
    }
    
    const std::array<int, 3> & rsize() const{
      return nr_;
    }

    const std::array<int, 3> & gsize() const{
      return ng_;
    }

    const math::d3vector & rspacing() const{
      return rspacing_;
    }

    const math::d3vector & gspacing() const{
      return gspacing_;
    }
    
    const math::d3vector & rlength() const{
      return rlength_;
    }

    const math::d3vector & glength() const{
      return glength_;
    }

    int rtotalsize() const {
      return nr_[0]*nr_[1]*nr_[2];
    }

    int gtotalsize() const {
      return ng_[0]*ng_[1]*ng_[2];
    }

		math::d3vector rvector(const int ix, const int iy, const int iz) const {
			return math::d3vector(ix*rspacing_[0], iy*rspacing_[1], iz*rspacing_[2]) - 0.5*rlength();
		}
		
		math::d3vector gvector(const int ix, const int iy, const int iz) const {
			math::d3vector g{ix*gspacing()[0], iy*gspacing()[1], iz*gspacing()[2]};
			for(int idir = 0; idir < 3; idir++) {
				if(g[idir] > 0.5*glength()[idir]) g[idir] -= glength()[idir];
			}
			return g;
		}

		bool g_is_zero(const int ix, const int iy, const int iz) const {
			return (ix == 0 and iy == 0 and iz == 0);
		}

		double g2(const int ix, const int iy, const int iz) const {
			return norm(gvector(ix, iy, iz));
		}

	private:
    double ecut_;
    std::array<int, 3> nr_;
    std::array<int, 3> ng_;

    math::d3vector rspacing_;
    math::d3vector gspacing_;
    
    math::d3vector rlength_;
    math::d3vector glength_;
    
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::plane_wave", "[basis]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
  {
    
    SECTION("Cubic cell"){

      ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

      double ecut = 20.0;
      
      basis::plane_wave pw(cell, ecut);

      REQUIRE(pw.ecut() == Approx(ecut));

      REQUIRE(pw.rtotalsize() == 8000);
      REQUIRE(pw.gtotalsize() == 8000);
      
      REQUIRE(pw.rspacing()[0] == 0.5_a);
      REQUIRE(pw.rspacing()[1] == 0.5_a);
      REQUIRE(pw.rspacing()[2] == 0.5_a);
      
      REQUIRE(pw.gspacing()[0] == 0.3141592654_a);
      REQUIRE(pw.gspacing()[1] == 0.3141592654_a);
      REQUIRE(pw.gspacing()[2] == 0.3141592654_a);
      
      REQUIRE(pw.rsize()[0] == 20);
      REQUIRE(pw.rsize()[1] == 20);
      REQUIRE(pw.rsize()[2] == 20);

      REQUIRE(pw.gsize()[0] == 20);
      REQUIRE(pw.gsize()[1] == 20);
      REQUIRE(pw.gsize()[2] == 20);

    }

    SECTION("Parallelepipedic cell"){

      ions::UnitCell cell(d3vector(77.7, 0.0, 0.0), d3vector(0.0, 14.14, 0.0), d3vector(0.0, 0.0, 23.25));

      double ecut = 37.9423091;
      
      basis::plane_wave pw(cell, ecut);

      REQUIRE(pw.ecut() == Approx(ecut));

      REQUIRE(pw.rtotalsize() == 536640);
      REQUIRE(pw.gtotalsize() == 536640);
	    
      REQUIRE(pw.rspacing()[0] == 0.3613953488_a);
      REQUIRE(pw.rspacing()[1] == 0.3625641026_a);
      REQUIRE(pw.rspacing()[2] == 0.36328125_a);
      
      REQUIRE(pw.gspacing()[0] == 0.0404323379_a);
      REQUIRE(pw.gspacing()[1] == 0.2221776983_a);
      REQUIRE(pw.gspacing()[2] == 0.1351222647_a);
      
      REQUIRE(pw.rsize()[0] == 215);
      REQUIRE(pw.rsize()[1] == 39);
      REQUIRE(pw.rsize()[2] == 64);

      REQUIRE(pw.gsize()[0] == 215);
      REQUIRE(pw.gsize()[1] == 39);
      REQUIRE(pw.gsize()[2] == 64);

    }

  }
}
#endif

    
#endif

