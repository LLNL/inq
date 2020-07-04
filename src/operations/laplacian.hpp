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

#include <multi/adaptors/fftw.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/cufft.hpp>
#endif

#include <cassert>

namespace inq {
namespace operations {
	
void laplacian_add(basis::field_set<basis::fourier_space, complex> const & ff, basis::field_set<basis::fourier_space, complex> const & laplff){
		
	//DATAOPERATIONS LOOP + GPU::RUN 4D
#ifdef ENABLE_CUDA
	
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						laplffcub = begin(laplff.cubic()),
						ffcub = begin(ff.cubic())]
					 __device__ (auto ist, auto iz, auto iy, auto ix){
						 
						 double lapl = -0.5*(-point_op.g2(ix, iy, iz));
						 laplffcub[ix][iy][iz][ist] += lapl*ffcub[ix][iy][iz][ist];
						 
					 });
	
#else
			auto point_op = ff.basis().point_op();
			
			for(int ix = 0; ix < laplff.basis().local_sizes()[0]; ix++){
				for(int iy = 0; iy < laplff.basis().local_sizes()[1]; iy++){
					for(int iz = 0; iz < laplff.basis().local_sizes()[2]; iz++){
						double lapl = -0.5*(-point_op.g2(ix, iy, iz));
						for(int ist = 0; ist < laplff.set_part().local_size(); ist++) laplff.cubic()[ix][iy][iz][ist] += lapl*ff.cubic()[ix][iy][iz][ist];
					}
				}
			}
			
#endif

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	void laplacian_in_place(basis::field_set<basis::fourier_space, complex> const & ff){


		//DATAOPERATIONS LOOP + GPU::RUN 4D
#ifdef ENABLE_CUDA

		gpu::run(ff.set_part().local_size(), ff.basis().local_sizes()[2], ff.basis().local_sizes()[1], ff.basis().local_sizes()[0],
						 [point_op = ff.basis().point_op(),
							ffcub = begin(ff.cubic())]
						 __device__ (auto ist, auto iz, auto iy, auto ix){
								 
							 double lapl = -0.5*(-point_op.g2(ix, iy, iz));
							 ffcub[ix][iy][iz][ist] = ffcub[ix][iy][iz][ist]*lapl;
								 
						 });

#else

		auto point_op = ff.basis().point_op();
		
		for(int ix = 0; ix < ff.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < ff.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < ff.basis().local_sizes()[2]; iz++){
					double lapl = -0.5*(-point_op.g2(ix, iy, iz));
					for(int ist = 0; ist < ff.set_part().local_size(); ist++) ff.cubic()[ix][iy][iz][ist] *= lapl;
				}
			}
		}
			
#endif
		
	}
	
}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("function operations::laplacian", "[operations::laplacian]") {

	using namespace inq;
	using namespace Catch::literals;

	
	
}


#endif

#endif

