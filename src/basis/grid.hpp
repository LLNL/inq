/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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

  class grid {

  public:

		const static int dimension = 3;
		
		template<class lattice_vectors_type>
		grid(const lattice_vectors_type & lattice_vectors, std::array<int, 3> nr) : nr_{nr}{
			for(int idir = 0; idir < 3; idir++){
				rlength_[idir] = length(lattice_vectors[idir]);
				ng_[idir] = nr_[idir];
				rspacing_[idir] = rlength_[idir]/nr_[idir];
				glength_[idir] = 2.0*M_PI/rspacing_[idir];
				gspacing_[idir] = glength_[idir]/ng_[idir];
			}

			npoints_ = nr_[0]*long(nr_[1])*nr_[2];
			
		}

    const std::array<int, 3> & rsize() const{
      return nr_;
    }

    const math::d3vector & rspacing() const{
      return rspacing_;
    }

		double diagonal_length() const {
			return length(rlength_);
		}
		
    const math::d3vector & rlength() const{
      return rlength_;
    }

    long size() const {
      return npoints_;
    }

		long num_points() const {
			return npoints_;
		}

		double volume_element() const {
			return rspacing_[0]*rspacing_[1]*rspacing_[2];
		}

		template <class output_stream>
    void info(output_stream & out) const {
      out << "PLANE WAVE BASIS SET:" << std::endl;
      out << "  Grid size   = " << rsize()[0] << " x " << rsize()[1] << " x " << rsize()[2] << std::endl;
			out << "  Spacing [b] = " << rspacing() << std::endl;
			out << std::endl;
    }

		friend auto sizes(const grid & gr){
			return gr.nr_;
		}
		
	protected:
		
    std::array<int, 3> nr_;
    std::array<int, 3> ng_;

    math::d3vector rspacing_;
    math::d3vector gspacing_;
    
    math::d3vector rlength_;
    math::d3vector glength_;

		long npoints_;
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::grid", "[grid]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
}
#endif

    
#endif
