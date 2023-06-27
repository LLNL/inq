/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__FOURIER_SPACE
#define INQ__BASIS__FOURIER_SPACE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include "grid.hpp"
#include <cassert>
#include <array>

namespace inq {
namespace basis {

class real_space;

  class fourier_space : public grid{

  public:

		using reciprocal_space = real_space;
		
    fourier_space(const grid & grid_basis):
			grid(grid_basis){
			
			cubic_part_ = {inq::parallel::partition(nr_[0]), inq::parallel::partition(nr_[1]), inq::parallel::partition(nr_[2], comm())};

			base::part_ = cubic_part_[2];
			base::part_ *= nr_[0]*long(nr_[1]);
			
			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_part_[idir].local_size();			
    }

		bool spherical() const {
			return spherical_g_grid_;
		}

		auto volume_element() const {
			return cell().volume()/(size()*size());
		}

		class point_operator {

		public:

			point_operator(std::array<int, 3> const & ng, vector3<double, covariant> const & gspacing, std::array<inq::parallel::partition, 3> const & dist, systems::cell::cell_metric metric):
				ng_(ng),
				gspacing_(gspacing),
				cubic_part_(dist),
				metric_(metric)
			{
			}

			GPU_FUNCTION auto gvector(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				
				//FFTW generates a grid from 0 to 2pi/h, so we convert it to a
				//grid from -pi/h to pi/h
				
				auto ii = grid::to_symmetric_range(ng_, ix, iy, iz);
				return vector3<double, covariant>{ii[0]*gspacing_[0], ii[1]*gspacing_[1], ii[2]*gspacing_[2]};
			}

			GPU_FUNCTION auto gvector(int ix, int iy, int iz) const {
				auto ixg = cubic_part_[0].local_to_global(ix);
				auto iyg = cubic_part_[1].local_to_global(iy);
				auto izg = cubic_part_[2].local_to_global(iz);
				
				return gvector(ixg, iyg, izg);
			}

			GPU_FUNCTION auto gvector_cartesian(int ix, int iy, int iz) const {
				return metric_.to_cartesian(gvector(ix, iy, iz));
			}
			
			GPU_FUNCTION auto outside_sphere(int ix, int iy, int iz) const {
				auto gvec = gvector(ix, iy, iz);
				gvec[0] /= ng_[0]/2.0;
				gvec[1] /= ng_[1]/2.0;
				gvec[2] /= ng_[2]/2.0;
				return metric_.norm(gvec) > 1.0;
			}
			
			GPU_FUNCTION const auto & gspacing() const{
				return gspacing_;
			}

			GPU_FUNCTION bool g_is_zero(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				return (ix.value() == 0 and iy.value() == 0 and iz.value() == 0);
			}

			GPU_FUNCTION bool g_is_zero(int ix, int iy, int iz) const {
				auto ixg = cubic_part_[0].local_to_global(ix);
				auto iyg = cubic_part_[1].local_to_global(iy);
				auto izg = cubic_part_[2].local_to_global(iz);
				
				return g_is_zero(ixg, iyg, izg);
			}
			
			GPU_FUNCTION double g2(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				return metric_.norm(gvector(ix, iy, iz));
			}
			
			GPU_FUNCTION double g2(int ix, int iy, int iz) const {
				return metric_.norm(gvector(ix, iy, iz));
			}

			GPU_FUNCTION auto & metric() const {
				return metric_;
			}

			GPU_FUNCTION auto to_symmetric_range(int ix, int iy, int iz) const {
				return grid::to_symmetric_range(ng_, ix, iy, iz);
			}
			
			GPU_FUNCTION auto from_symmetric_range(vector3<int> ii) const {
				return grid::from_symmetric_range(ng_, ii);
			}
			
		private:
			
			std::array<int, 3> ng_;
			vector3<double, covariant> gspacing_;
			std::array<inq::parallel::partition, 3> cubic_part_;
			systems::cell::cell_metric metric_;
			
		};

		auto point_op() const {
			return point_operator{ng_, covspacing_, cubic_part_, cell_.metric()};			
		}

		template <typename ReciprocalBasis = reciprocal_space>
		auto reciprocal() const {
			return ReciprocalBasis(*this);
		}
		
	private:
		
  };

}
}
#endif

#ifdef INQ_BASIS_FOURIER_SPACE_UNIT_TEST
#undef INQ_BASIS_FOURIER_SPACE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

}
#endif
