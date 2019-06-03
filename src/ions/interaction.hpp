#ifndef IONS_INTERACTION
#define IONS_INTERACTION

#include <config.h>

#include <math/d3vector.hpp>
#include <cmath>

namespace ions {

  template <class array_cell, class array_charge, class array_positions, class array_forces>
  void interaction_energy(const int natoms, const double alpha,
			  const array_cell & cell, const array_charge & charge, const array_positions & positions, 
			  double & energy, array_forces & forces){

    using math::d3vector;

    energy = 0.0;
    for(int iatom = 0; iatom < natoms; iatom++) forces[iatom] = d3vector(0.0, 0.0, 0.0);
    
    for(int iatom = 0; iatom < natoms; iatom++){

      d3vector xi = positions[iatom];
      double zi = charge[iatom];
      
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
}

#ifdef UNIT_TEST
#include <catch.hpp>

#include <vector>
#include <valarray>

TEST_CASE("Function ions::interaction_energy", "[interaction_energy]") {

  using namespace Catch::literals;
  using math::d3vector;

  double aa = 7.653;
  
  std::vector<d3vector> cell(3);
  cell[0] = d3vector(aa, 0.0, 0.0);
  cell[1] = d3vector(0.0, aa, 0.0);
  cell[2] = d3vector(0.0, 0.0, aa);

  std::valarray<double> charge(4);
  charge = 3.0;

  std::vector<d3vector> positions(3);
  positions[0] = d3vector(0.0,    0.0,    0.0);
  positions[1] = d3vector(aa/2.0, aa/2.0, 0.0);
  positions[2] = d3vector(aa/2.0, 0.0,    aa/2.0);
  positions[3] = d3vector(0.0,    aa/2.0, aa/2.0);

  double energy;
  std::vector<d3vector> forces(4);

  ions::interaction_energy(4, 0.21, cell, charge, positions, energy, forces);

  REQUIRE(energy == -9.99517178_a);
  
}
#endif

#endif

// Local Variables:
// mode: c++
// End:
