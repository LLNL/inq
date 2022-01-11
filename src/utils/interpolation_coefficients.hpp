/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__INTERPOLATION_COEFFICIENTS
#define INQ__UTILS__INTERPOLATION_COEFFICIENTS

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

namespace inq {
namespace utils {


template <class PointsType>
auto interpolation_coefficients(PointsType const & points, double const pos = 0.0){

  int nn = points.size();

  PointsType coeff(nn);

  for(auto ii = 0; ii < nn; ii++){
    coeff[ii] = 1.0;
    for(auto kk = 0; kk < nn; kk++){
      if(ii == kk) continue;
      coeff[ii] = coeff[ii]*(pos - points[kk])/(points[ii] - points[kk]);
    }
  }

  return coeff;
}

}
}

#ifdef INQ_UTILS_INTERPOLATION_COEFFICIENTS_UNIT_TEST
#undef INQ_UTILS_INTERPOLATION_COEFFICIENTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("function utils::interpolation_coefficients", "utils") {

	using namespace inq;
	using namespace Catch::literals;

  SECTION("first order"){
    
    auto coeff = utils::interpolation_coefficients(std::vector<double>{-1.0, 1.0});
    
    CHECK(coeff[0] == 0.5_a);
    CHECK(coeff[1] == 0.5_a);
    CHECK(coeff[0] + coeff[1] == 1.0_a);
  }

  SECTION("first order shifted"){
    
    auto coeff = utils::interpolation_coefficients(std::vector<double>{0.0, 1.0}, 0.5);
    
    CHECK(coeff[0] == 0.5_a);
    CHECK(coeff[1] == 0.5_a);
    CHECK(coeff[0] + coeff[1] == 1.0_a);    
  }

  SECTION("third order"){
    
    auto coeff = utils::interpolation_coefficients(std::vector<double>{-2.0, -1.0, 1.0, 2.0});
    
    CHECK(coeff[0] == -0.166666666667_a);
    CHECK(coeff[1] ==  0.666666666667_a);
    CHECK(coeff[2] ==  0.666666666667_a);
    CHECK(coeff[3] == -0.166666666667_a);
    CHECK(coeff[0] + coeff[1] + coeff[2] + coeff[3] == 1.0_a);
  }
  
}

#endif

#endif


