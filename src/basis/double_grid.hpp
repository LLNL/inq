/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__DOUBLE_GRID
#define INQ__BASIS__DOUBLE_GRID

/*
  Copyright (C) 2021 Xavier Andrade

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

#include <cassert>

#include <mpi3/environment.hpp>

#include <math/array.hpp>
#include <math/vector3.hpp>
#include <utils/interpolation_coefficients.hpp>
#include <utils/partition.hpp>

namespace inq {
namespace basis {

class double_grid {
	
public:


  template <class Function>
  auto value(Function const & func, math::vector3<double> pos){
    return func(pos);    
  }
  
  
};


}
}


#ifdef INQ_BASIS_DOUBLE_GRID_UNIT_TEST
#undef INQ_BASIS_DOUBLE_GRID_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::double_grid", "[basis::double_grid]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;


  basis::double_grid dg;

  CHECK(dg.value([](auto point){ return 1.0; }, {1.0, 2.0, 3.0}) == 1.0);

  
  
}
#endif

#endif
