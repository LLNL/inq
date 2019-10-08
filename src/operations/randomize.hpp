/* -*- indent-tabs-mode: t -*- */

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

#include <basis/field_set.hpp>
#include <cstdlib>

namespace operations {

	template <class field_set_type>
  void randomize(field_set_type & phi){
		srand48(0);

		//DATAOPERATIONS LOOP 2D (can be flattened, random number generator)
		for(int kk = 0; kk < phi.basis().size(); kk++) {
			for(int ii = 0; ii < phi.set_size(); ii++)	phi[ii][kk] = complex(drand48(), 0.0);
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
