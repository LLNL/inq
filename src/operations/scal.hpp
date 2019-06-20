/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__SCAL
#define OPERATIONS__SCAL

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

  template <class array_1d>
  void scal_invsqrt(const states::ks_states st, const basis::plane_wave & basis, const array_1d & factor, states::coefficients & phi){
    
    assert(size(factor) == st.num_states());
    
    for(int kk = 0; kk < basis.num_points(); kk++) {
      for(int ii = 0; ii < st.num_states(); ii++) phi.linear[ii][kk] /= sqrt(factor[ii]);
    }
    
  }
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::scal", "[scal]") {

	using namespace Catch::literals;
}


#endif

#endif
