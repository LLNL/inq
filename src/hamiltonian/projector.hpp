#ifndef HAMILTONIAN_PROJECTOR
#define HAMILTONIAN_PROJECTOR

#include <multi/array.hpp>

#include <math/d3vector.hpp>
#include <math/spherical_harmonic.hpp>
#include <pseudo/pseudopotential.hpp>
#include <ions/unitcell.hpp>
#include <basis/plane_wave.hpp>

namespace hamiltonian {

  class projector {

  public:
    projector(const basis::plane_wave & basis, const ions::UnitCell & cell, pseudo::pseudopotential ps, math::d3vector atom_position):
      sphere_(basis, cell, atom_position, ps.projector_radius()), nproj_(ps.num_projectors_lm()), matrix_({nproj_, sphere_.size()}) {

      std::vector<double> grid(sphere_.size()), proj(sphere_.size());

      // calculate the distance to the atom for each point
      for(int ipoint = 0; ipoint < sphere_.size(); ipoint++) grid[ipoint] = length(basis.rvector(sphere_.points()[ipoint]) - atom_position);

      for(int iproj = 0; iproj < ps.num_projectors_l(); iproj++){
	// interpolate the value of the radial part of the projectors to the sphere points
	ps.projector(iproj).value(sphere_.size(), grid, proj);

	int l = ps.projector_l(iproj);

	// now construct the projector with the spherical harmonics
	for(int m = -l; l <= m; m++){
	  for(int ipoint = 0; ipoint < sphere_.size(); ipoint++){
	    auto point = basis.rvector(sphere_.points()[ipoint]) - atom_position;
	    matrix_[iproj][ipoint] = proj[ipoint]*math::spherical_harmonic(l, m, point);
	  }
	}
	
      }

      
    }

    void apply(){

    }
    
    
  private:

    basis::spherical_grid sphere_;
    int nproj_;
    boost::multi::array<double, 2> matrix_;
    
  };
  
}

#endif
