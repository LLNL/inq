/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__FOURIER_SPACE
#define INQ__BASIS__FOURIER_SPACE

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

#include <inq_config.h>

#include <math/vec3d.hpp>
#include <ions/unitcell.hpp>
#include "grid.hpp"
#include <cassert>
#include <array>

namespace inq {
namespace basis {

  class fourier_space : public grid{

  public:
		
    fourier_space(const grid & grid_basis):
			grid(grid_basis){
			
			cubic_dist_ = {inq::utils::partition(nr_[0]), inq::utils::partition(nr_[1]), inq::utils::partition(nr_[2], comm())};

			base::part_ = cubic_dist_[2];
			base::part_ *= nr_[0]*long(nr_[1]);
			
			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_dist_[idir].local_size();			
    }

		bool spherical() const {
			return spherical_g_grid_;
		}

		auto volume_element() const {
			return rspacing_[0]*rspacing_[1]*rspacing_[2]/size();
		}

		class point_operator {

		public:

			point_operator(std::array<int, 3> const & ng, math::vec3d const & gspacing, math::vec3d const & glength):
				ng_(ng),
				gspacing_(gspacing),
				glength_(glength){
			}

			GPU_FUNCTION math::vec3d gvector(const int ix, const int iy, const int iz) const {
				
				//FFTW generates a grid from 0 to 2pi/h, so we convert it to a
				//grid from -pi/h to pi/h
				
				auto ii = to_symmetric_range(ng_, ix, iy, iz);
				return math::vec3d{ii[0]*gspacing_[0], ii[1]*gspacing_[1], ii[2]*gspacing_[2]};
			}
			
			GPU_FUNCTION const math::vec3d & glength() const{
				return glength_;
			}
			
			GPU_FUNCTION auto radius() const {
				return 0.5*std::min({glength_[0], glength_[1], glength_[2]});
			}
			
			GPU_FUNCTION auto outside_sphere(const double g2) const {
				return g2 > radius()*radius();
			}		
			
			GPU_FUNCTION const math::vec3d & gspacing() const{
				return gspacing_;
			}
			
			bool g_is_zero(const int ix, const int iy, const int iz) const {
				return (ix == 0 and iy == 0 and iz == 0);
			}
			
			GPU_FUNCTION double g2(const int ix, const int iy, const int iz) const {
				return norm(gvector(ix, iy, iz));
			}
			
		private:
			
			std::array<int, 3> ng_;
			math::vec3d gspacing_;
			math::vec3d glength_;
			
		};

		auto point_op() const {
			return point_operator(ng_, gspacing_, glength_);
		}
		
	private:
		
  };

}
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::fourier_space", "[fourier_space]") {

}
#endif

    
#endif
