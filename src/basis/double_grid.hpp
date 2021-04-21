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

  double_grid(bool enabled, int order = 5):
    enabled_(enabled),
    min_(-order + 1),
    max_(order),
    npoints_(max_ - min_ + 1)
  {

    math::array<double, 1> points(npoints_);
    for(auto ii = 0; ii < npoints_; ii++) points[ii] = min_ + ii;
    coeff_ = utils::interpolation_coefficients(points, 1.0/3.0);

  }

  template <class Function>
  auto value(Function const & func, math::vector3<double> pos, math::vector3<double> spacing) const {

    if(not enabled_) return func(pos);

    decltype(func(pos)) val = 0.0;
    
    for(int k0 = min_; k0 <= max_; k0++){
      for(int k1 = min_; k1 <= max_; k1++){
        for(int k2 = min_; k2 <= max_; k2++){
          
          double fac = coeff_[k0 - min_]*coeff_[k1 - min_]*coeff_[k2 - min_]/27.0;
          
          for(int i0 = -1; i0 <= 1; i0++){
            for(int i1 = -1; i1 <= 1; i1++){
              for(int i2 = -1; i2 <= 1; i2++){
                  
                val += fac*func(pos + spacing*math::vector3<double>{i0*(1.0 - 3.0*k0), i1*(1.0 - 3.0*k1), i2*(1.0 - 3.0*k2)});
                  
              }
            }
          }
          
          
        }
      }
    }
    
    return val;
  }

  auto spacing_factor() const {
    if(not enabled_) return 1.0;
    return 3.0;
  }

  auto & enabled() const {
    return enabled_;
  }
  
private:

  bool enabled_;
  int min_;
  int max_;
  int npoints_;
  math::array<double, 1> coeff_;
  
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


  basis::double_grid dg({0.3, 0.3, 0.3});

  CHECK(dg.value([](auto point){ return 1.0; }, {1.0, 2.0, 3.0}) == 1.0_a);

  
  
}
#endif

#endif
