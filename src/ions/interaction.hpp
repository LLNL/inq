#ifndef IONS_INTERACTION
#define IONS_INTERACTION

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
#include <limits>
#include <complex>

namespace ions {

  template <class cell_type, class array_charge, class array_positions, class array_forces>
  void interaction_energy(const int natoms, const double alpha,
			  const cell_type & cell, const array_charge & charge, const array_positions & positions, 
			  double & energy, array_forces & forces){

    using math::d3vector;

    energy = 0.0;
    for(int iatom = 0; iatom < natoms; iatom++) forces[iatom] = d3vector(0.0, 0.0, 0.0);

    double rcut = 6.0/alpha;

    for(int iatom = 0; iatom < natoms; iatom++){
      double zi = charge[iatom];
      
      periodic_replicas rep(cell, positions[iatom], rcut);

      for(int irep = 0; irep < rep.size(); irep++){
	d3vector xi = rep[irep];
	
	for(int jatom = 0; jatom < natoms; jatom++){
	  double zj = charge[jatom];
	  
	  d3vector rij = xi - positions[jatom];
	  double rr = sqrt(norm(rij));
	  
	  if(rr < 1.0e-5) continue;
	  
	  double eor = erfc(alpha*rr)/rr;
	  
	  energy += 0.5*zi*zj*eor;
	  forces[jatom] -= zi*zj*rij*(eor + 2.0*alpha/sqrt(M_PI)*exp(-pow(alpha*rr, 2))/(rr*rr));
	}
      }
      
    }

    // self-interaction
    double total_charge = 0.0;
    for(int iatom = 0; iatom < natoms; iatom++){
      double zi = charge[iatom];

      total_charge += zi;
      energy -= alpha*zi*zi/sqrt(M_PI);
    }

    // G = 0 energy
    energy -= M_PI*total_charge*total_charge/(2.0*alpha*alpha*cell.volume());

    double gcut = std::numeric_limits<double>::max();
    for(int idir = 0; idir < 3; idir++) std::min(gcut, norm(cell.b(idir)));
    gcut = sqrt(gcut);
      
    const int isph = ceil(9.5*alpha/gcut);

    std::vector<std::complex<double> > phase(natoms);
    
    for(int ix = -isph; ix <= isph; ix++){
      for(int iy = -isph; iy <= isph; iy++){
	for(int iz = -isph; iz <= isph; iz++){
	  
	  const int ss = ix*ix + iy*iy + iz*iz;

	  if(ss == 0 || ss > isph*isph) continue;

	  d3vector gg = ix*cell.b(0) + iy*cell.b(1) + iz*cell.b(2);
	  double gg2 = norm(gg);

	  double exparg = -0.25*gg2/(alpha*alpha);

	  if(exparg < -36.0) continue;

	  double factor = 2.0*M_PI/cell.volume()*exp(exparg)/gg2;

	  std::complex<double> sumatoms = 0.0;
	  for(int iatom = 0; iatom < natoms; iatom++){
	    double gx = gg*positions[iatom];
	    auto aa = charge[iatom]*std::complex<double>(cos(gx), sin(gx));
	    phase[iatom] = aa;
	    sumatoms += aa;
	  }

	  energy += factor*std::real(sumatoms*std::conj(sumatoms));
	  
	  for(int iatom = 0; iatom < natoms; iatom++){
	    for(int idir = 0; idir < 3; idir++){
	      std::complex<double> tmp = std::complex<double>(0.0, 1.0)*gg[idir]*phase[iatom];
	      forces[iatom][idir] -= factor*std::real(std::conj(tmp)*sumatoms + tmp*std::conj(sumatoms));
	    }
	  }
	  
	}
      }
    }

    //forces are not properly validated right now
    for(int iatom = 0; iatom < natoms; iatom++) forces[iatom] = d3vector(0.0, 0.0, 0.0);
    
  }
}

#ifdef UNIT_TEST
#include <catch.hpp>

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
    
    double energy;
    std::vector<d3vector> forces(4);
    
    ions::interaction_energy(4, 0.21, cell, charge, positions, energy, forces);
    
    REQUIRE(energy == -10.78368187_a); //this number comes from Octopus
    
  }

  SECTION("Diamond"){

    double aa = 6.74065308785213;
    
    ions::UnitCell cell(d3vector(0.0, aa/2.0, aa/2.0), d3vector(aa/2.0, 0.0, aa/2.0), d3vector(aa/2.0, aa/2.0, 0.0));

    const double charge[2] = {4.0, 4.0};
    
    std::vector<d3vector> positions(2);
    positions[0] = cell.crystal_to_cart(d3vector(0.0,  0.0,  0.0 ));
    positions[1] = cell.crystal_to_cart(d3vector(0.25, 0.25, 0.25));
    
    double energy;
    std::vector<d3vector> forces(2);

    ions::interaction_energy(2, 0.21, cell, charge, positions, energy, forces);

    REQUIRE(energy == -12.78641217_a); //this number comes from Octopus

  }

  SECTION("Iron"){
    
    double aa = 5.3970578;
    
    ions::UnitCell cell(d3vector(-aa/2.0, aa/2.0, aa/2.0), d3vector(aa/2.0, -aa/2.0, aa/2.0), d3vector(aa/2.0, aa/2.0, -aa/2.0));

    const double charge = 16.0;
    
    const d3vector position(0.0, 0.0, 0.0);
    
    double energy;
    std::vector<d3vector> forces(1);

    ions::interaction_energy(1, 0.21, cell, &charge, &position, energy, forces);

    REQUIRE(energy == -86.31033718_a); //this number comes from Octopus
    
  }
    
}
#endif

#endif

// Local Variables:
// mode: c++
// End:
