/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__LAPLACIAN
#define OPERATIONS__LAPLACIAN

/*
 Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa.

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

#include <gpu/run.hpp>
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <operations/space.hpp>

#include <multi/adaptors/fftw.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/cufft.hpp>
#endif

#include <cassert>

namespace inq {
namespace operations {
	
void laplacian_add(basis::field_set<basis::fourier_space, complex> const & ff, basis::field_set<basis::fourier_space, complex> const & laplff){
		
	//DATAOPERATIONS GPU::RUN 4D
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						laplffcub = begin(laplff.cubic()),
						ffcub = begin(ff.cubic()),
						cubic_dist_0 = ff.basis().cubic_dist(0),
						cubic_dist_1 = ff.basis().cubic_dist(1),
						cubic_dist_2 = ff.basis().cubic_dist(2)
						]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){

						 auto ixg = cubic_dist_0.local_to_global(ix);
						 auto iyg = cubic_dist_1.local_to_global(iy);
						 auto izg = cubic_dist_2.local_to_global(iz);
						 
						 double lapl = -0.5*(-point_op.g2(ixg, iyg, izg));
						 laplffcub[ix][iy][iz][ist] += lapl*ffcub[ix][iy][iz][ist];
						 
					 });
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void laplacian_in_place(basis::field_set<basis::fourier_space, complex> const & ff){
	
	//DATAOPERATIONS GPU::RUN 4D
	gpu::run(ff.set_part().local_size(), ff.basis().local_sizes()[2], ff.basis().local_sizes()[1], ff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						ffcub = begin(ff.cubic()),
						cubic_dist_0 = ff.basis().cubic_dist(0),
						cubic_dist_1 = ff.basis().cubic_dist(1),
						cubic_dist_2 = ff.basis().cubic_dist(2)
						]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){

						 auto ixg = cubic_dist_0.local_to_global(ix);
						 auto iyg = cubic_dist_1.local_to_global(iy);
						 auto izg = cubic_dist_2.local_to_global(iz);
						 
						 double lapl = -0.5*(-point_op.g2(ixg, iyg, izg));
						 ffcub[ix][iy][iz][ist] = ffcub[ix][iy][iz][ist]*lapl;
						 
					 });
}

}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST
#undef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("function operations::laplacian", "[operations::laplacian]") {

	using namespace inq;
	using namespace Catch::literals;

	
	
}


#endif

#endif

