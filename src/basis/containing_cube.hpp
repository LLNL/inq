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
void containing_cube(const BasisType & parent_grid, PosType const & pos, double radius, math::vector3<int> & lo, math::vector3<int> & hi){
	for(int idir = 0; idir < 3; idir++){
		
		//this doesnt work for non-orthogonal cells yet
		if(not parent_grid.cell().is_cartesian()){
			lo[idir] = parent_grid.symmetric_range_begin(idir);
			hi[idir] = parent_grid.symmetric_range_end(idir);
			continue;
		}
		
		lo[idir] = floor((pos[idir] - radius)/parent_grid.rspacing()[idir]) - 1;
		hi[idir] = ceil((pos[idir] + radius)/parent_grid.rspacing()[idir]) + 1;
		
		lo[idir] = std::clamp(lo[idir], parent_grid.symmetric_range_begin(idir), parent_grid.symmetric_range_end(idir));
		hi[idir] = std::clamp(hi[idir], parent_grid.symmetric_range_begin(idir), parent_grid.symmetric_range_end(idir));
	}
}

}
}

#ifdef INQ_BASIS_CONTAINING_CUBE_UNIT_TEST
#undef INQ_BASIS_CONTAINING_CUBE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("class basis::containing_cube", "[basis::containing_cube]") {
	
}
#endif

#endif

