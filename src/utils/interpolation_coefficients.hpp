/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__INTERPOLATION_COEFFICIENTS
#define INQ__UTILS__INTERPOLATION_COEFFICIENTS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
#endif

#ifdef INQ_UTILS_INTERPOLATION_COEFFICIENTS_UNIT_TEST
#undef INQ_UTILS_INTERPOLATION_COEFFICIENTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

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


