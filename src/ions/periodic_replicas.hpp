/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__PERIODIC_REPLICAS
#define INQ__IONS__PERIODIC_REPLICAS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/vector3.hpp>

#include <utils/profiling.hpp>

#include <vector>
#include <cmath>

namespace inq {
namespace ions {

class periodic_replicas{

public:

	template <class cell_array>
	periodic_replicas(const cell_array & cell, vector3<double> position, const double range){

		CALI_CXX_MARK_FUNCTION;
		
		assert(range >= 0.0);
		
		position = cell.position_in_cell(position);

		replicas_.push_back(position);		
		
		vector3<int> neigh_max{0, 0, 0};
		//we should use floor here, but since we check later, round is more reliable
		for(int idir = 0; idir < cell.periodicity(); idir++) neigh_max[idir] = round(range/sqrt(norm(cell[idir])));
		
		for(int ix = -neigh_max[0]; ix <= neigh_max[0]; ix++){
			for(int iy = -neigh_max[1]; iy <= neigh_max[1]; iy++){
				for(int iz = -neigh_max[2]; iz <= neigh_max[2]; iz++){

					if(ix == 0 and iy == 0 and iz == 0) continue;
					
					vector3<double> reppos = position + ix*cell[0] + iy*cell[1] + iz*cell[2];
            
					if(norm(reppos - position) <= range*range + 0.1) replicas_.push_back(reppos);
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

	std::vector<vector3<double>> replicas_;

};
  
}
}
#endif

#ifdef INQ_IONS_PERIODIC_REPLICAS_UNIT_TEST
#undef INQ_IONS_PERIODIC_REPLICAS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	
  {
    systems::cell cell(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 10.0, 0.0), vector3<double>(0.0, 0.0, 10.0));

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

		SECTION("Periodicity 2"){

			systems::cell cell2(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 10.0, 0.0), vector3<double>(0.0, 0.0, 10.0), 2);

			CHECK(cell2.periodicity() == 2);
			
			ions::periodic_replicas rep(cell2, vector3<double>(5.0, 5.0, 5.0), 11.0);
      
      CHECK(rep.size() == 5);

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
      CHECK(rep[3][1] ==  5.0_a);
      CHECK(rep[3][2] == -5.0_a);

      CHECK(rep[4][0] ==  5.0_a);
      CHECK(rep[4][1] == -5.0_a);
      CHECK(rep[4][2] == -5.0_a);

    }

		SECTION("Periodicity 0"){

			systems::cell cell0(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 10.0, 0.0), vector3<double>(0.0, 0.0, 10.0), 0);

			CHECK(cell0.periodicity() == 0);
			
			ions::periodic_replicas rep(cell0, vector3<double>(5.0, 5.0, 5.0), 11.0);
      
      CHECK(rep.size() == 1);

      CHECK(rep[0][0] == -5.0_a);
      CHECK(rep[0][1] == -5.0_a);
      CHECK(rep[0][2] == -5.0_a);
    }
  }

}
#endif
