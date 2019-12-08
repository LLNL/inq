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

#include <pcg-cpp/pcg_random.hpp>

namespace operations {

	template <class field_set_type>
  void randomize(field_set_type & phi){

		auto seed = phi.basis().size()*phi.set_size();
		
		pcg32 rng(seed);
		
		std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

		for(uint64_t ix = 0; ix < uint64_t(phi.basis().rsize()[0]); ix++){
			for(uint64_t iy = 0; iy < uint64_t(phi.basis().rsize()[1]); iy++){
				for(uint64_t iz = 0; iz < uint64_t(phi.basis().rsize()[2]); iz++){
					for(uint64_t ist = 0; ist < uint64_t(phi.dist().local_size()); ist++) {
						phi.cubic()[ix][iy][iz][ist] = uniform_dist(rng);
					}
				}
			}
		}
		
  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::randomize", "[operations::randomize]") {

	using namespace Catch::literals;

}


#endif

#endif
