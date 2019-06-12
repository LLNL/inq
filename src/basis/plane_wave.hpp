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
				glength_[idir] = 2.0*M_PI/rspacing_[idir];
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
				glength_[idir] = 2.0*M_PI/rspacing_[idir];
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

		double diagonal_length() const {
			return length(rlength_);
		}
		
    const math::d3vector & rlength() const{
      return rlength_;
    }

    const math::d3vector & glength() const{
      return glength_;
    }

    long size() const {
      return nr_[0]*long(nr_[1])*nr_[2];
    }
		
    long rtotalsize() const {
      return nr_[0]*long(nr_[1])*nr_[2];
    }

    long gtotalsize() const {
      return ng_[0]*long(ng_[1])*ng_[2];
    }

		math::d3vector rvector(const int ix, const int iy, const int iz) const {

			// shift so that the 0 0 0 point is at the center of the grid
			return math::d3vector(ix*rspacing_[0], iy*rspacing_[1], iz*rspacing_[2]) - 0.5*rlength();
		}

		template <class int_array>
		math::d3vector rvector(const int_array & indices) const {
			return rvector(indices[0], indices[1], indices[2]);
		}

		double r2(const int ix, const int iy, const int iz) const {
			return norm(rvector(ix, iy, iz));
		}
		
		math::d3vector gvector(const int ix, const int iy, const int iz) const {

			//FFTW generates a grid from 0 to 2pi/h, so we convert it to a
			//grid from -pi/h to pi/h
			
			math::d3vector g{ix*gspacing()[0], iy*gspacing()[1], iz*gspacing()[2]};
			for(int idir = 0; idir < 3; idir++) {
				if(g[idir] >= 0.5*glength()[idir]) g[idir] -= glength()[idir];
			}
			return g;
		}

		bool g_is_zero(const int ix, const int iy, const int iz) const {
			return (ix == 0 and iy == 0 and iz == 0);
		}

		double g2(const int ix, const int iy, const int iz) const {
			return norm(gvector(ix, iy, iz));
		}

		double volume_element() const {
			return rspacing_[0]*rspacing_[1]*rspacing_[2];
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
      
      REQUIRE(pw.gspacing()[0] == 0.6283185307_a);
      REQUIRE(pw.gspacing()[1] == 0.6283185307_a);
      REQUIRE(pw.gspacing()[2] == 0.6283185307_a);
      
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
      
      REQUIRE(pw.gspacing()[0] == 0.0808646758_a);
      REQUIRE(pw.gspacing()[1] == 0.4443553965_a);
      REQUIRE(pw.gspacing()[2] == 0.2702445293_a);
      
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

// Local variables:
// eval:(setq indent-tabs-mode t tab-width 2)
// End:
