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
		
    real_space(const ions::UnitCell & cell, const input::basis & basis_input, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			grid(cell, calculate_dimensions(cell, basis_input), basis_input.spherical_grid(), cell.periodic_dimensions(), comm){
    }

		real_space(const grid & grid_basis, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			grid(grid_basis){
			
			cubic_dist_ = {utils::partition(nr_[0], comm), utils::partition(nr_[1]), utils::partition(nr_[2])};

			base::part_ = cubic_dist_[0];
			base::part_ *= nr_[1]*long(nr_[2]);
			
			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_dist_[idir].local_size();		
    }

		GPU_FUNCTION math::vec3d rvector(const int ix, const int iy, const int iz) const {
			auto ii = this->to_symmetric_range(ix, iy, iz);
			return math::vec3d{ii[0]*rspacing()[0], ii[1]*rspacing()[1], ii[2]*rspacing()[2]};
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

		auto enlarge(int factor, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()) const {
			return real_space(grid(cell_.enlarge(factor), {factor*nr_[0], factor*nr_[1], factor*nr_[2]}, spherical_g_grid_, periodic_dimensions_, comm));
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

      CHECK(rs.size() == 8000);
      
      CHECK(rs.rspacing()[0] == 0.5_a);
      CHECK(rs.rspacing()[1] == 0.5_a);
      CHECK(rs.rspacing()[2] == 0.5_a);
      
      CHECK(rs.sizes()[0] == 20);
      CHECK(rs.sizes()[1] == 20);
      CHECK(rs.sizes()[2] == 20);

    }

    SECTION("Parallelepipedic cell"){

      ions::UnitCell cell(vec3d(77.7, 0.0, 0.0), vec3d(0.0, 14.14, 0.0), vec3d(0.0, 0.0, 23.25));

      double ecut = 37.9423091;
      
      basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

      CHECK(rs.size() == 536640);
	    
      CHECK(rs.rspacing()[0] == 0.3613953488_a);
      CHECK(rs.rspacing()[1] == 0.3625641026_a);
      CHECK(rs.rspacing()[2] == 0.36328125_a);
      
      CHECK(rs.sizes()[0] == 215);
      CHECK(rs.sizes()[1] == 39);
			CHECK(rs.sizes()[2] == 64);

			auto rs3x = rs.enlarge(3);
			
      CHECK(rs3x.rspacing()[0] == 0.3613953488_a);
      CHECK(rs3x.rspacing()[1] == 0.3625641026_a);
      CHECK(rs3x.rspacing()[2] == 0.36328125_a);
      
      CHECK(rs3x.sizes()[0] == 3*215);
      CHECK(rs3x.sizes()[1] == 3*39);
			CHECK(rs3x.sizes()[2] == 3*64);
			
    }

  }
}
#endif

    
#endif
