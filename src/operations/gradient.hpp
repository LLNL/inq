/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__GRADIENT
#define OPERATIONS__GRADIENT

/*
 Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa, Alexey Karstev.

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

#include <basis/field.hpp>
#include <cassert>

namespace operations {
	auto gradient(basis::field<basis::fourier_space, complex> const & ff){
		basis::field_set<basis::fourier_space, complex> grad(ff.basis(), 3);
		for(int ix = 0; ix < ff.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < ff.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < ff.basis().sizes()[2]; iz++){
					auto gvec = ff.basis().gvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) grad.cubic()[ix][iy][iz][idir] = complex(0.0, 1.0)*gvec[idir]*ff.cubic()[ix][iy][iz];
				}
			}
		}
		return grad;
	}
	auto gradient(basis::field<basis::real_space, complex> const & ff){
		auto ff_fourier = operations::space::to_fourier(ff); 			// To Fourier space
		auto grad_fourier = gradient(ff_fourier); 				// Computer gradient in Fourier space
		auto grad_real = operations::space::to_real(grad_fourier, false); 	// To real space
		return grad_real;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("function operations::gradient", "[operations::gradient]") {

//	using namespace Catch::literals;
// take sin function
}


#endif

#endif

