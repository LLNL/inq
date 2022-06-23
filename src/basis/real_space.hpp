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


#include <cassert>
#include <array>

#include <basis/grid.hpp>
#include <gpu/run.hpp>
#include <ions/unitcell.hpp>
#include <math/vector3.hpp>
#include <systems/box.hpp>

namespace inq {
namespace basis {

  class real_space : public grid {

  public:
		
    real_space(systems::box const & box, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			grid(box, calculate_dimensions(box), box.spherical_grid_value(), box.double_grid_value(), box.periodic_dimensions_value(), comm)
		{
    }

		real_space(const grid & grid_basis):
			grid(grid_basis){
			
			cubic_dist_ = {inq::parallel::partition(nr_[0], grid_basis.comm()), inq::parallel::partition(nr_[1]), inq::parallel::partition(nr_[2])};

			base::part_ = cubic_dist_[0];
			base::part_ *= nr_[1]*long(nr_[2]);
			
			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_dist_[idir].local_size();		
    }
		
		real_space(real_space && old, boost::mpi3::communicator & new_comm):
			real_space(grid(grid(old), new_comm))
		{
		}

		class point_operator {

		public:

			point_operator(std::array<int, 3> const & nr, math::vector3<double> const & rspacing, std::array<inq::parallel::partition, 3> const & dist, ions::UnitCell::cell_metric metric):
				nr_(nr),
				rspacing_(rspacing),
				cubic_dist_(dist),
				metric_(metric)
			{
			}

			GPU_FUNCTION auto to_symmetric_range(int ix, int iy, int iz) const {
				return grid::to_symmetric_range(nr_, ix, iy, iz);
			}
			
			GPU_FUNCTION auto from_symmetric_range(math::vector3<int> ii) const {
				return grid::from_symmetric_range(nr_, ii);
			}

			GPU_FUNCTION auto rvector(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				auto ii = grid::to_symmetric_range(nr_, ix, iy, iz);
				return math::vector3<double, math::contravariant>{ii[0]*rspacing_[0], ii[1]*rspacing_[1], ii[2]*rspacing_[2]};
			}
			
			GPU_FUNCTION auto rvector(int ix, int iy, int iz) const {
				auto ixg = cubic_dist_[0].local_to_global(ix);
				auto iyg = cubic_dist_[1].local_to_global(iy);
				auto izg = cubic_dist_[2].local_to_global(iz);
				
				return rvector(ixg, iyg, izg);
			}

			GPU_FUNCTION auto rvector_cartesian(int ix, int iy, int iz) const {
				return metric_.to_cartesian(rvector(ix, iy, iz));
			}

			GPU_FUNCTION auto rvector_cartesian(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				return metric_.to_cartesian(rvector(ix, iy, iz));
			}
			
			template <class int_array>
			GPU_FUNCTION auto rvector(const int_array & indices) const {
				return rvector(indices[0], indices[1], indices[2]);
			}
			
			template <typename IndexType>
			GPU_FUNCTION double r2(IndexType ix, IndexType iy, IndexType iz) const {
				return metric_.norm(rvector(ix, iy, iz));
			}

			template <typename IndexType>
			GPU_FUNCTION double rlength(IndexType ix, IndexType iy, IndexType iz) const {
				return metric_.length(rvector(ix, iy, iz));
			}
			
			GPU_FUNCTION auto & cubic_dist() const {
				return cubic_dist_;
			}

			GPU_FUNCTION auto & metric() const {
				return metric_;
			}				
			
		private:
			
			std::array<int, 3> nr_;
			math::vector3<double> rspacing_;
			std::array<inq::parallel::partition, 3> cubic_dist_;
			ions::UnitCell::cell_metric metric_;
			
		};
			
		friend auto operator==(const real_space & rs1, const real_space & rs2){
			bool equal = rs1.nr_[0] == rs2.nr_[0] and rs1.nr_[1] == rs2.nr_[1] and rs1.nr_[2] == rs2.nr_[2];
			equal = equal and rs1.rspacing()[0] == rs2.rspacing()[0];
			equal = equal and rs1.rspacing()[1] == rs2.rspacing()[1];
			equal = equal and rs1.rspacing()[2] == rs2.rspacing()[2];
			return equal;
		}

		auto enlarge(int factor) const {
			return real_space(grid(cell_.enlarge(factor), {factor*nr_[0], factor*nr_[1], factor*nr_[2]}, spherical_g_grid_, double_grid_.enabled(), periodic_dimensions_, this->comm()));
		}

		auto refine(double factor, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()) const {
			assert(factor > 0.0);
			return real_space(grid(cell_, {(int) round(factor*nr_[0]), (int) round(factor*nr_[1]), (int) round(factor*nr_[2])}, spherical_g_grid_,  double_grid_.enabled(), periodic_dimensions_, comm));
		}
		
		auto volume_element() const {
			return rspacing_[0]*rspacing_[1]*rspacing_[2];
		}

		auto gcutoff() const {
			auto max_spacing = std::max({rspacing_[0], rspacing_[1],rspacing_[2]});

			return double_grid_.spacing_factor()*M_PI/max_spacing;
		}

		auto point_op() const {
			return point_operator(nr_, rspacing_, cubic_dist_, cell_.metric());
		}
		
	private:

		static std::array<int, 3> calculate_dimensions(systems::box const & box){
			std::array<int, 3> nr;
			double spacing = box.spacing_value();
			
			// make the spacing conmensurate with the grid
			// OPTIMIZATION: we can select a good size here for the FFT
			for(int idir = 0; idir < 3; idir++){
				double rlength = length(box[idir]);
				nr[idir] = round(rlength/spacing);
			}
			
			return nr;
		}

  };

}
}

#ifdef INQ_BASIS_REAL_SPACE_UNIT_TEST
#undef INQ_BASIS_REAL_SPACE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::real_space", "[basis::real_space]") {

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

	using math::vector3;

	auto comm = boost::mpi3::environment::get_world_instance();
	
  {
    
    SECTION("Cubic cell"){

			systems::box box = systems::box::cubic(10.0_b).cutoff_energy(20.0_Ha);
      basis::real_space rs(box, comm);

      CHECK(rs.size() == 8000);
      
      CHECK(rs.rspacing()[0] == 0.5_a);
      CHECK(rs.rspacing()[1] == 0.5_a);
      CHECK(rs.rspacing()[2] == 0.5_a);
      
      CHECK(rs.sizes()[0] == 20);
      CHECK(rs.sizes()[1] == 20);
      CHECK(rs.sizes()[2] == 20);
			
			basis::real_space new_rs(basis::real_space(rs), boost::mpi3::environment::get_self_instance());
			
			CHECK(rs.sizes() == new_rs.sizes());
			CHECK(new_rs.local_sizes() == new_rs.sizes());	
			CHECK(rs.periodic_dimensions() == new_rs.periodic_dimensions());
			
    }

    SECTION("Parallelepipedic cell"){

			systems::box box = systems::box::orthorhombic(77.7_b, 14.14_b, 23.25_b).cutoff_energy(37.9423091_Ha);
      basis::real_space rs(box, comm);

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
