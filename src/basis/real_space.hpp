/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__REAL_SPACE
#define INQ__BASIS__REAL_SPACE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>
#include <array>

#include <basis/grid.hpp>
#include <gpu/run.hpp>
#include <math/vector3.hpp>

namespace inq {
namespace basis {

class fourier_space;

  class real_space : public grid {

  public:

		using reciprocal_space = fourier_space;

		real_space(systems::cell const & cell, double const & spacing, parallel::communicator comm):
			grid(cell, calculate_dimensions(cell, spacing), comm)
		{
    }

		real_space(const grid & grid_basis):
			grid(grid_basis){
			
			cubic_part_ = {inq::parallel::partition(nr_[0], grid_basis.comm()), inq::parallel::partition(nr_[1]), inq::parallel::partition(nr_[2])};

			base::part_ = cubic_part_[0];
			base::part_ *= nr_[1]*long(nr_[2]);
			
			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_part_[idir].local_size();		
    }
		
		real_space(real_space && old, parallel::communicator & new_comm):
			real_space(grid(grid(old), new_comm))
		{
		}

		real_space(real_space && old, parallel::communicator && new_comm):
			real_space(grid(grid(old), new_comm))
		{
		}

		class point_operator {

		public:

			point_operator(std::array<int, 3> const & nr, vector3<double, contravariant> const & rspacing, std::array<inq::parallel::partition, 3> const & dist, systems::cell::cell_metric metric):
				nr_(nr),
				rspacing_(rspacing),
				cubic_part_(dist),
				metric_(metric)
			{
			}

			GPU_FUNCTION auto to_symmetric_range(int ix, int iy, int iz) const {
				return grid::to_symmetric_range(nr_, ix, iy, iz);
			}
			
			GPU_FUNCTION auto from_symmetric_range(vector3<int> ii) const {
				return grid::from_symmetric_range(nr_, ii);
			}

			GPU_FUNCTION auto rvector(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				auto ii = grid::to_symmetric_range(nr_, ix, iy, iz);
				return vector3<int, contravariant>{ii[0], ii[1], ii[2]}*rspacing_;
			}
			
			GPU_FUNCTION auto rvector(int ix, int iy, int iz) const {
				auto ixg = cubic_part_[0].local_to_global(ix);
				auto iyg = cubic_part_[1].local_to_global(iy);
				auto izg = cubic_part_[2].local_to_global(iz);
				
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
			
			GPU_FUNCTION auto & cubic_part(int idim) const {
				return cubic_part_[idim];
			}
			
			GPU_FUNCTION auto & metric() const {
				return metric_;
			}				
			
		private:
			
			std::array<int, 3> nr_;
			vector3<double, contravariant> rspacing_;
			std::array<inq::parallel::partition, 3> cubic_part_;
			systems::cell::cell_metric metric_;
			
		};
			
		friend auto operator==(const real_space & rs1, const real_space & rs2){
			bool equal = rs1.nr_[0] == rs2.nr_[0] and rs1.nr_[1] == rs2.nr_[1] and rs1.nr_[2] == rs2.nr_[2];
			equal = equal and rs1.rspacing()[0] == rs2.rspacing()[0];
			equal = equal and rs1.rspacing()[1] == rs2.rspacing()[1];
			equal = equal and rs1.rspacing()[2] == rs2.rspacing()[2];
			return equal;
		}

		auto enlarge(int factor) const {
			return real_space(grid(cell_.enlarge(factor), {factor*nr_[0], factor*nr_[1], factor*nr_[2]}, this->comm()));
		}

		auto enlarge(vector3<int> factor) const {
			return real_space(grid(cell_.enlarge(factor), {factor[0]*nr_[0], factor[1]*nr_[1], factor[2]*nr_[2]},  this->comm()));
		}
		
		auto refine(double factor) const {
			assert(factor > 0.0);
			return real_space(grid(cell_, {(int) round(factor*nr_[0]), (int) round(factor*nr_[1]), (int) round(factor*nr_[2])}, this->comm()));
		}
		
		auto volume_element() const {
			return cell().volume()/size();
		}

		auto gcutoff() const {
			auto max_spacing = std::max({rspacing_[0], rspacing_[1],rspacing_[2]});

			return M_PI/max_spacing;
		}

		auto point_op() const {
			return point_operator(nr_, conspacing_, cubic_part_, cell_.metric());
		}

		template <typename ReciprocalBasis = reciprocal_space>
		auto reciprocal() const {
			return ReciprocalBasis(*this);
		}

		static auto gcutoff(systems::cell const & cell, double const & spacing){
			auto nr = calculate_dimensions(cell, spacing);

			auto max_spacing = 0.0;
			for(int idir = 0; idir < 3; idir++){
				auto actual_spacing = length(cell[idir])/nr[idir];
				max_spacing = std::max(max_spacing, actual_spacing);
			}

			return M_PI/max_spacing;
		}
		
	private:

		static std::array<int, 3> calculate_dimensions(systems::cell const & cell, double const & spacing){
			std::array<int, 3> nr;
			
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
#endif

#ifdef INQ_BASIS_REAL_SPACE_UNIT_TEST
#undef INQ_BASIS_REAL_SPACE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

	
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
  {
    
    SECTION("Cubic cell"){

      basis::real_space rs(systems::cell::cubic(10.0_b), /* spacing = */ 0.49672941, comm);

      CHECK(rs.size() == 8000);

			CHECK(rs.cell().volume() == rs.volume_element()*rs.size());
      
      CHECK(rs.rspacing()[0] == 0.5_a);
      CHECK(rs.rspacing()[1] == 0.5_a);
      CHECK(rs.rspacing()[2] == 0.5_a);
      
      CHECK(rs.sizes()[0] == 20);
      CHECK(rs.sizes()[1] == 20);
      CHECK(rs.sizes()[2] == 20);

			CHECK(rs.gcutoff() == 6.2831853072_a);
			CHECK(rs.gcutoff() == basis::real_space::gcutoff(systems::cell::cubic(10.0_b), /* spacing = */ 0.49672941));
			
			basis::real_space new_rs(basis::real_space(rs), parallel::communicator{boost::mpi3::environment::get_self_instance()});
			
			CHECK(rs.sizes() == new_rs.sizes());
			CHECK(new_rs.local_sizes() == new_rs.sizes());	
			CHECK(rs.cell().periodicity() == new_rs.cell().periodicity());

			CHECK(rs.part().rank() == comm.rank());
			CHECK(rs.cubic_part(0).rank() == comm.rank()%rs.cubic_part(0).comm_size());
			CHECK(rs.cubic_part(1).rank() == comm.rank()%rs.cubic_part(1).comm_size());
			CHECK(rs.cubic_part(2).rank() == comm.rank()%rs.cubic_part(2).comm_size());
			
			rs.shift();
			CHECK(rs.part().rank() == (comm.rank() + 1)%comm.size());
			CHECK(rs.cubic_part(0).rank() == (comm.rank() + 1)%rs.cubic_part(0).comm_size());
			CHECK(rs.cubic_part(1).rank() == (comm.rank() + 1)%rs.cubic_part(1).comm_size());
			CHECK(rs.cubic_part(2).rank() == (comm.rank() + 1)%rs.cubic_part(2).comm_size());			
			
    }

    SECTION("Parallelepipedic cell"){

      basis::real_space rs(systems::cell::orthorhombic(77.7_b, 14.14_b, 23.25_b), /*spacing =*/ 0.36063925, comm);

      CHECK(rs.size() == 536640);

			CHECK(rs.cell().volume() == rs.volume_element()*rs.size());
			
      CHECK(rs.rspacing()[0] == 0.3613953488_a);
      CHECK(rs.rspacing()[1] == 0.3625641026_a);
      CHECK(rs.rspacing()[2] == 0.36328125_a);
      
      CHECK(rs.sizes()[0] == 215);
      CHECK(rs.sizes()[1] == 39);
			CHECK(rs.sizes()[2] == 64);

			CHECK(rs.gcutoff() == 8.6478249389_a);
			CHECK(rs.gcutoff() == basis::real_space::gcutoff(systems::cell::orthorhombic(77.7_b, 14.14_b, 23.25_b), /*spacing =*/ 0.36063925));
			
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

    SECTION("Non-orthogonal cell"){

			auto a = 3.567095_A;
      basis::real_space rs(systems::cell::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}), /*spacing = */ a/30.0_b, comm);

      CHECK(rs.size() == 9261);

			CHECK(rs.cell().volume() == Approx(0.25*pow(a.in_atomic_units(), 3)));
			CHECK(rs.cell().volume() == rs.volume_element()*rs.size());
			
      CHECK(rs.rspacing()[0] == 0.2269756405_a);
      CHECK(rs.rspacing()[1] == 0.2269756405_a);
      CHECK(rs.rspacing()[2] == 0.2269756405_a);
      
      CHECK(rs.sizes()[0] == 21);
      CHECK(rs.sizes()[1] == 21);
			CHECK(rs.sizes()[2] == 21);

			CHECK(rs.gcutoff() == 13.8411005126_a);
			CHECK(rs.gcutoff() == basis::real_space::gcutoff(systems::cell::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}), /*spacing = */ a/30.0_b));

    }

  }
}
#endif
