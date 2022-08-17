/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__PERIODIC_REPLICAS
#define INQ__IONS__PERIODIC_REPLICAS

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

#include <math/vector3.hpp>

#include <utils/profiling.hpp>

#include <vector>
#include <cmath>

namespace inq {
namespace ions {

class periodic_replicas{

public:

	template <class cell_array>
	periodic_replicas(const cell_array & cell, math::vector3<double> position, const double range){

		CALI_CXX_MARK_FUNCTION;
		
		assert(range >= 0.0);
		
		position = cell.position_in_cell(position);

		replicas_.push_back(position);		
		
		math::vector3<int> neigh_max(3);
		//we should use floor here, but since we check later, round is more reliable
		for(int idir = 0; idir < 3; idir++) neigh_max[idir] = round(range/sqrt(norm(cell[idir]))); 
      
		for(int ix = -neigh_max[0]; ix <= neigh_max[0]; ix++){
			for(int iy = -neigh_max[1]; iy <= neigh_max[1]; iy++){
				for(int iz = -neigh_max[2]; iz <= neigh_max[2]; iz++){

					if(ix == 0 and iy == 0 and iz == 0) continue;
					
					math::vector3<double> reppos = position + ix*cell[0] + iy*cell[1] + iz*cell[2];
            
					if(norm(reppos - position) <= range*range) replicas_.push_back(reppos);
				}
			}
		}
      
	}

	const auto & operator[](const int i) const {
		return replicas_[i];
	}

	auto size(){
		return replicas_.size();
	}
    
private:

	std::vector<math::vector3<double>> replicas_;

};    
  
}
}


#ifdef INQ_IONS_PERIODIC_REPLICAS_UNIT_TEST
#undef INQ_IONS_PERIODIC_REPLICAS_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("class ions::periodic_replicas", "[periodic_replicas]") {

	using namespace inq;
	using namespace Catch::literals;
	using math::vector3;

  {
    ions::unit_cell cell(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 10.0, 0.0), vector3<double>(0.0, 0.0, 10.0));

    SECTION("Cubic cell 0"){
      ions::periodic_replicas rep(cell, vector3<double>(5.0, 5.0, 5.0), 9.5);
      
      CHECK(rep.size() == 1);
      
      CHECK(rep[0][0] == -5.0_a);
      CHECK(rep[0][1] == -5.0_a);
      CHECK(rep[0][2] == -5.0_a);

    }
    
    SECTION("Cubic cell 1"){
      ions::periodic_replicas rep(cell, vector3<double>(5.0, 5.0, 5.0), 10.0);
      
      CHECK(rep.size() == 7);

      CHECK(rep[0][0] == -5.0_a);
      CHECK(rep[0][1] == -5.0_a);
      CHECK(rep[0][2] == -5.0_a);
      
      CHECK(rep[1][0] == -15.0_a);
      CHECK(rep[1][1] == -5.0_a);
      CHECK(rep[1][2] == -5.0_a);
      
      CHECK(rep[2][0] == -5.0_a);
      CHECK(rep[2][1] == -15.0_a);
      CHECK(rep[2][2] == -5.0_a);
      
      CHECK(rep[3][0] == -5.0_a);
      CHECK(rep[3][1] == -5.0_a);
      CHECK(rep[3][2] == -15.0_a);
      
      CHECK(rep[4][0] == -5.0_a);
      CHECK(rep[4][1] == -5.0_a);
      CHECK(rep[4][2] == 5.0_a);
      
      CHECK(rep[5][0] == -5.0_a);
      CHECK(rep[5][1] == 5.0_a);
      CHECK(rep[5][2] == -5.0_a);

      CHECK(rep[6][0] == 5.0_a);
      CHECK(rep[6][1] == -5.0_a);
      CHECK(rep[6][2] == -5.0_a);

    }
    
    SECTION("Cubic cell 2"){
      ions::periodic_replicas rep(cell, vector3<double>(5.0, 5.0, 5.0), 11.0);
      
      CHECK(rep.size() == 7);

      CHECK(rep[0][0] == -5.0_a);
      CHECK(rep[0][1] == -5.0_a);
      CHECK(rep[0][2] == -5.0_a);
			
      CHECK(rep[1][0] == -15.0_a);
      CHECK(rep[1][1] == -5.0_a);
      CHECK(rep[1][2] == -5.0_a);
      
      CHECK(rep[2][0] == -5.0_a);
      CHECK(rep[2][1] == -15.0_a);
      CHECK(rep[2][2] == -5.0_a);
      
      CHECK(rep[3][0] == -5.0_a);
      CHECK(rep[3][1] == -5.0_a);
      CHECK(rep[3][2] == -15.0_a);
      
      CHECK(rep[4][0] == -5.0_a);
      CHECK(rep[4][1] == -5.0_a);
      CHECK(rep[4][2] == 5.0_a);
      
      CHECK(rep[5][0] == -5.0_a);
      CHECK(rep[5][1] == 5.0_a);
      CHECK(rep[5][2] == -5.0_a);

      CHECK(rep[6][0] == 5.0_a);
      CHECK(rep[6][1] == -5.0_a);
      CHECK(rep[6][2] == -5.0_a);

    }

    SECTION("Cubic cell 3"){
      ions::periodic_replicas rep(cell, vector3<double>(5.0, 5.0, 5.0), 15.0);
      
      CHECK(rep.size() == 19);
    }

    SECTION("Cubic cell 4"){
      ions::periodic_replicas rep(cell, vector3<double>(5.0, 5.0, 5.0), 18.0);
      
      CHECK(rep.size() == 27);
    }
    
    SECTION("Cubic cell 5"){
			
      ions::periodic_replicas rep(cell, vector3<double>(35.0, -205.0, 2035.0), 10.0);
      
      CHECK(rep.size() == 7);

      CHECK(rep[0][0] == -5.0_a);
      CHECK(rep[0][1] == -5.0_a);
      CHECK(rep[0][2] == -5.0_a);
			
      CHECK(rep[1][0] == -15.0_a);
      CHECK(rep[1][1] == -5.0_a);
      CHECK(rep[1][2] == -5.0_a);
      
      CHECK(rep[2][0] == -5.0_a);
      CHECK(rep[2][1] == -15.0_a);
      CHECK(rep[2][2] == -5.0_a);
      
      CHECK(rep[3][0] == -5.0_a);
      CHECK(rep[3][1] == -5.0_a);
      CHECK(rep[3][2] == -15.0_a);
      
      CHECK(rep[4][0] == -5.0_a);
      CHECK(rep[4][1] == -5.0_a);
      CHECK(rep[4][2] == 5.0_a);
      
      CHECK(rep[5][0] == -5.0_a);
      CHECK(rep[5][1] == 5.0_a);
      CHECK(rep[5][2] == -5.0_a);

      CHECK(rep[6][0] == 5.0_a);
      CHECK(rep[6][1] == -5.0_a);
      CHECK(rep[6][2] == -5.0_a);

    }
		
  }

}


#endif


#endif
