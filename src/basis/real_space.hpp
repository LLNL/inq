/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__REAL_SPACE
#define INQ__BASIS__REAL_SPACE

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
#include <ions/unitcell.hpp>
#include "grid.hpp"
#include <cassert>
#include <array>
#include <input/basis.hpp>
#include <gpu/run.hpp>

namespace inq {
namespace basis {

  class real_space : public grid {

  public:
		
    real_space(const ions::UnitCell & cell, const input::basis & basis_input, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			grid(cell, calculate_dimensions(cell, basis_input), basis_input.spherical_grid(), cell.periodic_dimensions(), comm){
    }

		real_space(const grid & grid_basis):
			grid(grid_basis){
			
			cubic_dist_ = {inq::utils::partition(nr_[0], grid_basis.comm()), inq::utils::partition(nr_[1]), inq::utils::partition(nr_[2])};

			base::part_ = cubic_dist_[0];
			base::part_ *= nr_[1]*long(nr_[2]);
			
			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_dist_[idir].local_size();		
    }

		GPU_FUNCTION math::vector3<double> rvector(utils::global_index ix, utils::global_index iy, utils::global_index iz) const {
			auto ii = this->to_symmetric_range(ix, iy, iz);
			return math::vector3<double>{ii[0]*rspacing()[0], ii[1]*rspacing()[1], ii[2]*rspacing()[2]};
		}

		GPU_FUNCTION math::vector3<double> rvector(int ix, int iy, int iz) const {
			auto ixg = cubic_dist(0).local_to_global(ix);
			auto iyg = cubic_dist(1).local_to_global(iy);
			auto izg = cubic_dist(2).local_to_global(iz);

			return rvector(ixg, iyg, izg);
		}
		
		template <class int_array>
		GPU_FUNCTION math::vector3<double> rvector(const int_array & indices) const {
			return rvector(indices[0], indices[1], indices[2]);
		}

		template <typename IndexType>
		double r2(IndexType ix, IndexType iy, IndexType iz) const {
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
			return real_space(grid(cell_.enlarge(factor), {factor*nr_[0], factor*nr_[1], factor*nr_[2]}, spherical_g_grid_, periodic_dimensions_, this->comm()));
		}

		auto refine(double factor, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()) const {
			assert(factor > 0.0);
			return real_space(grid(cell_, {(int) round(factor*nr_[0]), (int) round(factor*nr_[1]), (int) round(factor*nr_[2])}, spherical_g_grid_, periodic_dimensions_, comm));
		}
		
		auto volume_element() const {
			return rspacing_[0]*rspacing_[1]*rspacing_[2];
		}

		auto gcutoff() const {
			auto max_spacing = std::max({rspacing_[0], rspacing_[1],rspacing_[2]});

			return M_PI/max_spacing;
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
}

#ifdef INQ_BASIS_REAL_SPACE_UNIT_TEST
#undef INQ_BASIS_REAL_SPACE_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::real_space", "[basis::real_space]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;

  {
    
    SECTION("Cubic cell"){

      ions::UnitCell cell(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 10.0, 0.0), vector3<double>(0.0, 0.0, 10.0));

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

      ions::UnitCell cell(vector3<double>(77.7, 0.0, 0.0), vector3<double>(0.0, 14.14, 0.0), vector3<double>(0.0, 0.0, 23.25));

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

			auto rs_fine = rs.refine(2);
			
      CHECK(rs_fine.rspacing()[0] == Approx(0.5*0.3613953488));
      CHECK(rs_fine.rspacing()[1] == Approx(0.5*0.3625641026));
      CHECK(rs_fine.rspacing()[2] == Approx(0.5*0.36328125));
      
      CHECK(rs_fine.sizes()[0] == 2*215);
      CHECK(rs_fine.sizes()[1] == 2*39);
			CHECK(rs_fine.sizes()[2] == 2*64);

			CHECK(rs == rs.refine(1));
			
			auto rs_155 = rs.refine(1.55);
			
      CHECK(rs_155.rspacing()[0] == Approx(0.2333333333));
      CHECK(rs_155.rspacing()[1] == Approx(0.2356666667));
      CHECK(rs_155.rspacing()[2] == Approx(0.2348484848));
      
      CHECK(rs_155.sizes()[0] == 333);
      CHECK(rs_155.sizes()[1] == 60);
			CHECK(rs_155.sizes()[2] == 99);

    }

  }
}
#endif

    
#endif
