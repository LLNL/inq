////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// LDAFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LDAFUNCTIONAL_HPP
#define LDAFUNCTIONAL_HPP

#include <vector>
#include <cassert>

namespace functionals {
  class LDAFunctional {

  public:
    static void xc_unpolarized(const double rh, double &ee, double &vv);
    static void xc_polarized(const double rh, double &ee, double &vv);
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class functionals::LDAFunctional", "[LDAFunctional]") {

  using namespace Catch::literals;

    std::vector<double> nn({1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0, 1e1, 1e2, 1e3});
    std::vector<double> ex(nn.size()), vx(nn.size());

    //these values are from the same routine, to avoid changes in the
    //results. But they haven't been validated to other source yet.
    
    SECTION("spin unpolarized"){
      
      for(unsigned ii = 0; ii < nn.size(); ii++){
	functionals::LDAFunctional::xc_unpolarized(nn[ii], ex[ii], vx[ii]);
      }
      
      REQUIRE(ex[0] == -0.0121328183_a);
      REQUIRE(vx[0] == -0.0159054085_a);
      
      REQUIRE(ex[1] == -0.0246687773_a);
      REQUIRE(vx[1] == -0.0322046187_a);
      
      REQUIRE(ex[2] == -0.0495735123_a);
      REQUIRE(vx[2] == -0.0644958518_a);
      
      REQUIRE(ex[3] == -0.0988616346_a);
      REQUIRE(vx[3] == -0.1284302277_a);
      
      REQUIRE(ex[4] == -0.1970983191_a);
      REQUIRE(vx[4] == -0.2564000609_a);
      
      REQUIRE(ex[5] == -0.3962482024_a);
      REQUIRE(vx[5] == -0.5175699500_a);
      
      REQUIRE(ex[6] == -0.8092221329_a);
      REQUIRE(vx[6] == -1.0635879059_a);
      
      REQUIRE(ex[7] == -1.6819692861_a);
      REQUIRE(vx[7] == -2.2215642308_a);
      
      REQUIRE(ex[8] == -3.5407734574_a);
      REQUIRE(vx[8] == -4.6932262077_a);
      
      REQUIRE(ex[9] == -7.5211172181_a);
      REQUIRE(vx[9] == -9.9930315852_a);
    }
    
    SECTION("spin polarized"){
      for(unsigned ii = 0; ii < nn.size(); ii++){
	functionals::LDAFunctional::xc_polarized(nn[ii], ex[ii], vx[ii]);
      }

      REQUIRE(ex[0] == -0.0122936534_a);
      REQUIRE(vx[0] == -0.0161617994_a);
      
      REQUIRE(ex[1] == -0.0253096199_a);
      REQUIRE(vx[1] == -0.0332259764_a);
      
      REQUIRE(ex[2] == -0.0519716782_a);
      REQUIRE(vx[2] == -0.0682116367_a);
      
      REQUIRE(ex[3] == -0.1068678070_a);
      REQUIRE(vx[3] == -0.1404217239_a);
      
      REQUIRE(ex[4] == -0.2209158888_a);
      REQUIRE(vx[4] == -0.2909428124_a);
      
      REQUIRE(ex[5] == -0.4603409307_a);
      REQUIRE(vx[5] == -0.6080094243_a);
      
      REQUIRE(ex[6] == -0.968035_a);
      REQUIRE(vx[6] == -1.28248_a);
      
      REQUIRE(ex[7] == -2.05265_a);
      REQUIRE(vx[7] == -2.72561_a);
      
      REQUIRE(ex[8] == -4.37814_a);
      REQUIRE(vx[8] == -5.82279_a);
      
      REQUIRE(ex[9] == -9.37581_a);
      REQUIRE(vx[9] == -12.4826_a);
    }

}

#endif

#endif

// Local Variables:
// mode: c++
// End:
