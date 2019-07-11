/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__SHIFT
#define OPERATIONS__SHIFT

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
#include <cassert>

namespace operations {

	template <class array_1d, class coefficients_set_type>
  void shift(const array_1d & factor, const coefficients_set_type & shift, coefficients_set_type & phi){
    
    assert(size(factor) == phi.set_size());

    //DATAOPERATIONS
    
    for(int ipoint = 0; ipoint < phi.basis().num_points(); ipoint++) {
      for(int ist = 0; ist < phi.set_size(); ist++) phi.linear[ipoint][ist] += factor[ist]*shift.linear[ipoint][ist];
    }
    
  }
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::shift", "[shift]") {

	using namespace Catch::literals;
}


#endif

#endif
