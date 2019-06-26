/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef BASIS_REAL_SPACE
#define BASIS_REAL_SPACE

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

#include <math/d3vector.hpp>
#include <ions/unitcell.hpp>
#include "grid.hpp"
#include <cassert>
#include <array>

namespace basis {

  class real_space : public grid{

  public:
		
		real_space(ions::UnitCell & cell, std::array<int, 3> nr):
			grid(cell, nr){
		}
		
    real_space(ions::UnitCell & cell, const double & ecut):
			grid(cell, ecut){
    }
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::real_space", "[real_space]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
  {
    
    SECTION("Cubic cell"){

      ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

      double ecut = 20.0;
      
      basis::grid pw(cell, ecut);

      REQUIRE(pw.ecut() == Approx(ecut));

      REQUIRE(pw.rtotalsize() == 8000);
      REQUIRE(pw.gtotalsize() == 8000);
      
      REQUIRE(pw.rspacing()[0] == 0.5_a);
      REQUIRE(pw.rspacing()[1] == 0.5_a);
      REQUIRE(pw.rspacing()[2] == 0.5_a);
      
      REQUIRE(pw.gspacing()[0] == 0.6283185307_a);
      REQUIRE(pw.gspacing()[1] == 0.6283185307_a);
      REQUIRE(pw.gspacing()[2] == 0.6283185307_a);
      
      REQUIRE(pw.rsize()[0] == 20);
      REQUIRE(pw.rsize()[1] == 20);
      REQUIRE(pw.rsize()[2] == 20);

      REQUIRE(pw.gsize()[0] == 20);
      REQUIRE(pw.gsize()[1] == 20);
      REQUIRE(pw.gsize()[2] == 20);

    }

    SECTION("Parallelepipedic cell"){

      ions::UnitCell cell(d3vector(77.7, 0.0, 0.0), d3vector(0.0, 14.14, 0.0), d3vector(0.0, 0.0, 23.25));

      double ecut = 37.9423091;
      
      basis::grid pw(cell, ecut);

      REQUIRE(pw.ecut() == Approx(ecut));

      REQUIRE(pw.rtotalsize() == 536640);
      REQUIRE(pw.gtotalsize() == 536640);
	    
      REQUIRE(pw.rspacing()[0] == 0.3613953488_a);
      REQUIRE(pw.rspacing()[1] == 0.3625641026_a);
      REQUIRE(pw.rspacing()[2] == 0.36328125_a);
      
      REQUIRE(pw.gspacing()[0] == 0.0808646758_a);
      REQUIRE(pw.gspacing()[1] == 0.4443553965_a);
      REQUIRE(pw.gspacing()[2] == 0.2702445293_a);
      
      REQUIRE(pw.rsize()[0] == 215);
      REQUIRE(pw.rsize()[1] == 39);
      REQUIRE(pw.rsize()[2] == 64);

      REQUIRE(pw.gsize()[0] == 215);
      REQUIRE(pw.gsize()[1] == 39);
      REQUIRE(pw.gsize()[2] == 64);

    }

  }
}
#endif

    
#endif
