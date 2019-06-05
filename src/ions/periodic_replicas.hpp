#ifndef IONS_PERIODIC_REPLICAS
#define IONS_PERIODIC_REPLICAS

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
#include <vector>
#include <cmath>

namespace ions {
  class periodic_replicas{

  public:

    enum class error {
      NEGATIVE_RANGE
    };
    template <class cell_array>
    periodic_replicas(const cell_array & cell, const math::d3vector & position, const double range){

      using math::d3vector;
      
      if(range < 0) throw error::NEGATIVE_RANGE;

      std::vector<int> neigh_max(3);

      //we should use floor here, but since we check later, round is more reliable
      for(int idir = 0; idir < 3; idir++) neigh_max[idir] = round(range/sqrt(norm(cell[0]))); 
      
      for(int ix = -neigh_max[0]; ix <= neigh_max[0]; ix++){
	for(int iy = -neigh_max[1]; iy <= neigh_max[1]; iy++){
	  for(int iz = -neigh_max[2]; iz <= neigh_max[2]; iz++){
	    d3vector reppos = position + ix*cell[0] + iy*cell[1] + iz*cell[2];

	    if(norm(reppos - position) <= range*range) replicas_.push_back(reppos);
	  }
	}
      }
      
    }

    const math::d3vector & operator[](const int i) const {
      return replicas_[i];
    }

    const int size(){
      return replicas_.size();
    }
    
  private:

    std::vector<math::d3vector> replicas_;

  };    
  
}


#ifdef UNIT_TEST
#include <catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class ions::periodic_replicas", "[periodic_replicas]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
  {
    ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

    SECTION("Negative range"){
      REQUIRE_THROWS(ions::periodic_replicas(cell, d3vector(0.0, 0.0, 0.0), -1.0));
    }
    
    SECTION("Cubic cell 0"){
      ions::periodic_replicas rep(cell, d3vector(5.0, 5.0, 5.0), 9.5);
      
      REQUIRE(rep.size() == 1);
      
      REQUIRE(rep[0][0] == 5.0_a);
      REQUIRE(rep[0][1] == 5.0_a);
      REQUIRE(rep[0][2] == 5.0_a);

    }
    
    SECTION("Cubic cell 1"){
      ions::periodic_replicas rep(cell, d3vector(5.0, 5.0, 5.0), 10.0);
      
      REQUIRE(rep.size() == 7);
      
      REQUIRE(rep[0][0] == -5.0_a);
      REQUIRE(rep[0][1] == 5.0_a);
      REQUIRE(rep[0][2] == 5.0_a);
      
      REQUIRE(rep[1][0] == 5.0_a);
      REQUIRE(rep[1][1] == -5.0_a);
      REQUIRE(rep[1][2] == 5.0_a);
      
      REQUIRE(rep[2][0] == 5.0_a);
      REQUIRE(rep[2][1] == 5.0_a);
      REQUIRE(rep[2][2] == -5.0_a);
      
      REQUIRE(rep[3][0] == 5.0_a);
      REQUIRE(rep[3][1] == 5.0_a);
      REQUIRE(rep[3][2] == 5.0_a);
      
      REQUIRE(rep[4][0] == 5.0_a);
      REQUIRE(rep[4][1] == 5.0_a);
      REQUIRE(rep[4][2] == 15.0_a);
      
      REQUIRE(rep[5][0] == 5.0_a);
      REQUIRE(rep[5][1] == 15.0_a);
      REQUIRE(rep[5][2] == 5.0_a);

      REQUIRE(rep[6][0] == 15.0_a);
      REQUIRE(rep[6][1] == 5.0_a);
      REQUIRE(rep[6][2] == 5.0_a);

    }
    
    SECTION("Cubic cell 2"){
      ions::periodic_replicas rep(cell, d3vector(5.0, 5.0, 5.0), 11.0);
      
      REQUIRE(rep.size() == 7);
      
      REQUIRE(rep[0][0] == -5.0_a);
      REQUIRE(rep[0][1] == 5.0_a);
      REQUIRE(rep[0][2] == 5.0_a);
      
      REQUIRE(rep[1][0] == 5.0_a);
      REQUIRE(rep[1][1] == -5.0_a);
      REQUIRE(rep[1][2] == 5.0_a);
      
      REQUIRE(rep[2][0] == 5.0_a);
      REQUIRE(rep[2][1] == 5.0_a);
      REQUIRE(rep[2][2] == -5.0_a);
      
      REQUIRE(rep[3][0] == 5.0_a);
      REQUIRE(rep[3][1] == 5.0_a);
      REQUIRE(rep[3][2] == 5.0_a);
      
      REQUIRE(rep[4][0] == 5.0_a);
      REQUIRE(rep[4][1] == 5.0_a);
      REQUIRE(rep[4][2] == 15.0_a);
      
      REQUIRE(rep[5][0] == 5.0_a);
      REQUIRE(rep[5][1] == 15.0_a);
      REQUIRE(rep[5][2] == 5.0_a);

      REQUIRE(rep[6][0] == 15.0_a);
      REQUIRE(rep[6][1] == 5.0_a);
      REQUIRE(rep[6][2] == 5.0_a);

    }

    SECTION("Cubic cell 3"){
      ions::periodic_replicas rep(cell, d3vector(5.0, 5.0, 5.0), 15.0);
      
      REQUIRE(rep.size() == 19);
    }

    SECTION("Cubic cell 4"){
      ions::periodic_replicas rep(cell, d3vector(5.0, 5.0, 5.0), 18.0);
      
      REQUIRE(rep.size() == 27);
    }

  }

}


#endif


#endif

// Local Variables:
// mode: c++
// End:
