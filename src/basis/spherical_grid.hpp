#ifndef SPHERICAL_GRID_HPP
#define SPHERICAL_GRID_HPP

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
#include <cassert>
#include <array>

namespace basis {

  class spherical_grid {

  public:
    template <class basis>
    spherical_grid(const basis & parent_grid, const ions::UnitCell & cell, const math::d3vector & center_point, const double radius){
      ions::periodic_replicas rep(cell, center_point, radius);

      
      for(int ix = 0; ix < parent_grid.rsize()[0]; ix++){
	for(int iy = 0; iy < parent_grid.rsize()[1]; iy++){
	  for(int iz = 0; iz < parent_grid.rsize()[2]; iz++){
	    auto rpoint = parent_grid.rvector(ix, iy, iz);

	    for(int irep = 0; rep.size(); irep++){
	      if(norm(rpoint - rep[irep]) <= radius*radius){ 
		points_.push_back({ix, iy, iz});
		break;
	      }
	    }
	    
	  }
	}
      }

    }

    long size() const {
      return points_.size();
    }
    
  private:
    
    std::vector<std::array<int, 3> > points_;
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::spherical_grid", "[spherical_grid]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
  {
    
    SECTION("Cubic cell"){

      ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

      double ecut = 20.0;
      
      basis::plane_wave pw(cell, ecut);

      basis::spherical_grid sphere(pw, cell, {0.0, 0.0, 0.0}, 5.0);
      

    }

    SECTION("Parallelepipedic cell"){

      ions::UnitCell cell(d3vector(77.7, 0.0, 0.0), d3vector(0.0, 14.14, 0.0), d3vector(0.0, 0.0, 23.25));

      double ecut = 37.9423091;
      
      basis::plane_wave pw(cell, ecut);


    }

  }
}
#endif

    
#endif

