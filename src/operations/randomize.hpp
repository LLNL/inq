/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__RANDOMIZE
#define OPERATIONS__RANDOMIZE

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

#include <states/coefficients.hpp>
#include <cstdlib>

namespace operations {
	
  void randomize(const states::ks_states st, const basis::plane_wave & basis, states::coefficients & phi){
		srand48(0);

		for(int kk = 0; kk < basis.num_points(); kk++) {
			for(int ii = 0; ii < st.num_states(); ii++)	phi.linear[ii][kk] = complex(drand48(), drand48());
    }

  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::randomize", "[randomize]") {

	using namespace Catch::literals;

}


#endif

#endif
