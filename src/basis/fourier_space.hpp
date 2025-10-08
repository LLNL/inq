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

	systems::cell cell_;

public:

	using reciprocal_space = real_space;
		
	fourier_space(real_space const & rs):
		grid(rs.sizes(), rs.comm(), /*par_dim = */ 0),
		cell_(rs.cell())
	{
	}

	constexpr static auto covspacing(int) {
		return 2.0*M_PI;
	}
	
	auto & cell() const {
		return cell_;
	}

	auto volume_element() const {
		return cell().volume()/(size()*size());
	}
		
	friend auto operator==(const fourier_space & fs1, const fourier_space & fs2){
		return fs1.cell_ == fs2.cell_ and fs1.sizes_[0] == fs2.sizes_[0] and fs1.sizes_[1] == fs2.sizes_[1] and fs1.sizes_[2] == fs2.sizes_[2];
	}
		
	class point_operator {
			
		std::array<int, 3> sizes_;
		std::array<inq::parallel::partition, 3> cubic_part_;
		systems::cell::cell_metric metric_;

	public:
		
		point_operator(std::array<int, 3> const & nr, std::array<inq::parallel::partition, 3> const & dist, systems::cell::cell_metric metric):
			sizes_(nr),
			cubic_part_(dist),
			metric_(metric)
		{
		}

		GPU_FUNCTION auto gvector(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
				
			//FFTW generates a grid from 0 to 2pi/h, so we convert it to a
			//grid from -pi/h to pi/h
				
			auto ii = grid::to_symmetric_range(sizes_, ix, iy, iz);
			return vector3<double, covariant>{ii[0]*covspacing(0), ii[1]*covspacing(1), ii[2]*covspacing(2)};
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
			auto ivec = grid::to_symmetric_range(sizes_, cubic_part_[0].local_to_global(ix), cubic_part_[1].local_to_global(iy), cubic_part_[2].local_to_global(iz));
			double xx = double(ivec[0])/sizes_[0];
			double yy = double(ivec[1])/sizes_[1];
			double zz = double(ivec[2])/sizes_[2];
			return xx*xx + yy*yy + zz*zz > 0.25;
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
			return grid::to_symmetric_range(sizes_, ix, iy, iz);
		}
			
		GPU_FUNCTION auto from_symmetric_range(vector3<int> ii) const {
			return grid::from_symmetric_range(sizes_, ii);
		}
						
		GPU_FUNCTION auto & cubic_part(int idim) const {
			return cubic_part_[idim];
		}

		GPU_FUNCTION auto local_contains(vector3<int> const & ii) const {
			bool contains = true;
			for(int idir = 0; idir < 3; idir++){
				contains = contains and cubic_part_[idir].contains(ii[idir]);
			}
			return contains;
		}
			
	};

	auto point_op() const {
		return point_operator{sizes_, cubic_part_, cell_.metric()};
	}

	template <typename ReciprocalBasis = reciprocal_space>
	auto reciprocal() const {
		return ReciprocalBasis(cell_, sizes_, comm_);
	}

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
