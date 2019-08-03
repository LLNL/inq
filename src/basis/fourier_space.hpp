/* -*- indent-tabs-mode: t -*- */

#ifndef BASIS_FOURIER_SPACE
#define BASIS_FOURIER_SPACE

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
#include <ions/unitcell.hpp>
#include "grid.hpp"
#include <cassert>
#include <array>

namespace basis {

  class fourier_space : public grid{

  public:
		
    fourier_space(const grid & grid_basis):
			grid(grid_basis){
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

    const math::d3vector & glength() const{
      return glength_;
    }

		auto radius() const {
			return 0.5*std::min(glength_[0], std::min(glength_[1], glength_[2]));
		}

		auto outside_sphere(const double g2) const {
			if(not spherical_g_grid_) return false;
			return g2 > radius()*radius();
		}		

    const math::d3vector & gspacing() const{
      return gspacing_;
    }

    const std::array<int, 3> & gsize() const{
      return ng_;
    }
		
		bool g_is_zero(const int ix, const int iy, const int iz) const {
			return (ix == 0 and iy == 0 and iz == 0);
		}

		double g2(const int ix, const int iy, const int iz) const {
			return norm(gvector(ix, iy, iz));
		}

		bool spherical() const {
			return spherical_g_grid_;
		}
		
	private:
		
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::fourier_space", "[fourier_space]") {

}
#endif

    
#endif
