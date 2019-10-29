/* -*- indent-tabs-mode: t -*- */

#ifndef IONS__INTERACTION
#define IONS__INTERACTION

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

#include <config.h>
#include <math/d3vector.hpp>
#include <cmath>
#include <ions/periodic_replicas.hpp>
#include <math/array.hpp>
#include <limits>
#include <complex>

namespace ions {


	template <class cell_type, class geometry_type>
	void interaction_energy(const math::erf_range_separation & sep, const cell_type & cell, const geometry_type & geo, double & energy, double & eself){
		math::array<double, 1> charges(geo.num_atoms());

		for(int ii = 0; ii < geo.num_atoms(); ii++) charges[ii] = geo.atoms()[ii].charge();

		interaction_energy(geo.num_atoms(), cell, charges, geo.coordinates(), sep, energy, eself);
	}

	template <class cell_type, class array_charge, class array_positions>
  void interaction_energy(const int natoms, const cell_type & cell, const array_charge & charge, const array_positions & positions, 
													const math::erf_range_separation & sep, double & energy, double & eself){

    using math::d3vector;

    energy = 0.0;

		// the short range interaction
    double rcut = cell.diagonal_length() + sep.short_range_potential_radius();

    for(int iatom = 0; iatom < natoms; iatom++){
      auto zi = charge[iatom];
      
      periodic_replicas rep(cell, positions[iatom], rcut);

      for(unsigned irep = 0; irep < rep.size(); irep++){
				auto xi = rep[irep];
				
				for(int jatom = 0; jatom < natoms; jatom++){
					auto zj = charge[jatom];
					
					auto rij = xi - positions[jatom];
					auto rr = length(rij);
					
					if(rr < 1.0e-5) continue;
					
					energy += 0.5*zi*zj*sep.short_range_potential(rr);
				}
      }
      
    }

		eself = 0.0;

		//the self-interaction correction of the long range part
    for(int iatom = 0; iatom < natoms; iatom++){
      auto zi = charge[iatom];

			eself -= zi*zi*sep.self_interaction();
		}
			
	}
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

#include <vector>
#include <valarray>
#include <ions/unitcell.hpp>

TEST_CASE("Function ions::interaction_energy", "[interaction_energy]") {

  //Note: here we use Octopus' results for validation. However,
  //normally Octopus add a pseudopotential energy term to the ion-ion
  //energy that we include elsewhere, and that was removed for this
  //comparison.
  
  using namespace Catch::literals;
  using math::d3vector;

	const math::erf_range_separation sep(0.625);

  SECTION("Aluminum cubic cell"){
  
    double aa = 7.653;
    
    ions::UnitCell cell(d3vector(aa, 0.0, 0.0), d3vector(0.0, aa, 0.0), d3vector(0.0, 0.0, aa));
    
    std::valarray<double> charge(4);
    charge = 3.0;
    
    std::vector<d3vector> positions(4);
    positions[0] = d3vector(0.0,    0.0,    0.0);
    positions[1] = d3vector(aa/2.0, aa/2.0, 0.0);
    positions[2] = d3vector(aa/2.0, 0.0,    aa/2.0);
    positions[3] = d3vector(0.0,    aa/2.0, aa/2.0);
    
    double energy, eself;
    std::vector<d3vector> forces(4);
    
    ions::interaction_energy(4, cell, charge, positions, sep, energy, eself);
    
		//    REQUIRE(energy == -10.78368187_a); //this number comes from Octopus
    
  }

  SECTION("Diamond"){

    double aa = 6.74065308785213;
    
    ions::UnitCell cell(d3vector(0.0, aa/2.0, aa/2.0), d3vector(aa/2.0, 0.0, aa/2.0), d3vector(aa/2.0, aa/2.0, 0.0));

    const double charge[2] = {4.0, 4.0};
    
    std::vector<d3vector> positions(2);
    positions[0] = cell.crystal_to_cart(d3vector(0.0,  0.0,  0.0 ));
    positions[1] = cell.crystal_to_cart(d3vector(0.25, 0.25, 0.25));
    
    double energy, eself;
    std::vector<d3vector> forces(2);

    ions::interaction_energy(2, cell, charge, positions, sep, energy, eself);

		//    REQUIRE(energy == -12.78641217_a); //this number comes from Octopus

  }

  SECTION("Iron"){
    
    double aa = 5.3970578;
    
    ions::UnitCell cell(d3vector(-aa/2.0, aa/2.0, aa/2.0), d3vector(aa/2.0, -aa/2.0, aa/2.0), d3vector(aa/2.0, aa/2.0, -aa/2.0));

    const double charge = 16.0;
    
    const d3vector position(0.0, 0.0, 0.0);
    
    double energy, eself;
    std::vector<d3vector> forces(1);

    ions::interaction_energy(1, cell, &charge, &position, sep, energy, eself);

		//    REQUIRE(energy == -86.31033718_a); //this number comes from Octopus
    
  }
    
}
#endif

#endif
