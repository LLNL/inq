/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__CALCULATE_DENSITY
#define OPERATIONS__CALCULATE_DENSITY

/*
 Copyright (C) 2019 Xavier Andrade

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
#include <basis/field_set.hpp>
#include <math/complex.hpp>
#include <cstdlib>

namespace operations {

  template<class occupations_array_type>
  auto calculate_density(const occupations_array_type & occupations, const basis::field_set<basis::real_space, complex> & phi){
    basis::field<basis::real_space, double> density(phi.basis());

    //DATAOPERATIONS LOOP 2D
    for(int ipoint = 0; ipoint < phi.basis().size(); ipoint++){
      for(int ist = 0; ist < phi.set_size(); ist++) density[ipoint] += occupations[ist]*norm(phi[ist][ipoint]);
    }

    return density;
  }
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::randomize", "[randomize]") {

	using namespace Catch::literals;

}


#endif

#endif
