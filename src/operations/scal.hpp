/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SCAL
#define INQ__OPERATIONS__SCAL

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <basis/field_set.hpp>
#include <cassert>

namespace inq {
namespace operations {

//TODO: this should receive the operation as argument
template <class array_1d, class field_set_type>
void scal_invsqrt(const array_1d & factor, field_set_type & phi){
    
	assert(size(factor) == phi.set_size());

	//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifdef ENABLE_CUDA
	gpu::run(phi.set_size(), phi.basis().num_points(),
					 [factor, phimat = begin(phi.matrix()), fac = begin(factor)] __device__
					 (auto ist, auto ipoint){
						 phimat[ist][ipoint] /= sqrt(fac[ist]);
					 });
#else
	for(int kk = 0; kk < phi.basis().num_points(); kk++) {
		for(int ii = 0; ii < phi.set_size(); ii++) phi.matrix()[ii][kk] /= sqrt(factor[ii]);
	}
#endif
		
}
  
}
}

#ifdef INQ_OPERATIONS_SCAL_UNIT_TEST
#undef INQ_OPERATIONS_SCAL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("function operations::scal", "[scal]") {

	using namespace inq;
	using namespace Catch::literals;
}


#endif

#endif
