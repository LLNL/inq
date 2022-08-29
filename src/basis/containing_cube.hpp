/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__CONTAINING_CUBE
#define INQ__BASIS__CONTAINING_CUBE

/*
 Copyright (C) 2019-2021 Xavier Andrade

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

#include <math/vector3.hpp>

namespace inq {
namespace basis {

//returns the cube that contains the sphere, this makes the initialization O(1) instead of O(N)
template <class BasisType, typename PosType>
void containing_cube(const BasisType & grid, PosType const & pos, double radius, math::vector3<int> & lo, math::vector3<int> & hi){
	for(int idir = 0; idir < 3; idir++){
		
		//this doesnt work for non-orthogonal cells yet
		if(not grid.cell().is_cartesian()){
			lo[idir] = grid.symmetric_range_begin(idir);
			hi[idir] = grid.symmetric_range_end(idir);
			continue;
		}
		
		lo[idir] = floor((pos[idir] - radius)/grid.rspacing()[idir]) - 1;
		hi[idir] = ceil((pos[idir] + radius)/grid.rspacing()[idir]) + 2;
		
		lo[idir] = std::clamp(lo[idir], grid.symmetric_range_begin(idir), grid.symmetric_range_end(idir));
		hi[idir] = std::clamp(hi[idir], grid.symmetric_range_begin(idir), grid.symmetric_range_end(idir));
	}
}

}
}

#ifdef INQ_BASIS_CONTAINING_CUBE_UNIT_TEST
#undef INQ_BASIS_CONTAINING_CUBE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>
#include <math/array.hpp>
#include <basis/real_space.hpp>

TEST_CASE("class basis::containing_cube", "[basis::containing_cube]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;

	using math::vector3;

	auto comm = boost::mpi3::environment::get_world_instance();

	SECTION("Orthogonal box"){
		systems::box box = systems::box::orthorhombic(12.0_b, 14.0_b, 16.0_b).cutoff_energy(45.0_Ha);
		basis::real_space rs(box, comm);

		auto center = math::vector3{3.0, 2.0, 1.0};
		auto radius = 3.0;

		math::vector3<int> lo, hi;
		containing_cube(rs, center, radius, lo, hi);

		CHECK(lo[0] == -1);
		CHECK(lo[1] == -4);
		CHECK(lo[2] == -7);
		CHECK(hi[0] == 18);
		CHECK(hi[1] == 17);
		CHECK(hi[2] == 14);		
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto ii = rs.point_op().to_symmetric_range(ix, iy, iz);
					
					if(ii[0] >= lo[0] and ii[0] < hi[0] and
						 ii[1] >= lo[1] and ii[1] < hi[1] and
						 ii[2] >= lo[2] and ii[2] < hi[2]) continue;

					auto dist2 = norm(center - rs.point_op().rvector_cartesian(ii[0], iy, iz));
					CHECK(dist2 > radius*radius);
				}
			}
		}

		
	}

	
	
}
#endif

#endif

