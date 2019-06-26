/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef HAMILTONIAN_PROJECTOR
#define HAMILTONIAN_PROJECTOR

#include <multi/array.hpp>

#include <math/d3vector.hpp>
#include <math/spherical_harmonic.hpp>
#include <pseudo/pseudopotential.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/grid.hpp>
#include <basis/spherical_grid.hpp>

namespace hamiltonian {

  class projector {

  public:
    projector(const basis::grid & basis, const ions::UnitCell & cell, pseudo::pseudopotential ps, math::d3vector atom_position):
      sphere_(basis, cell, atom_position, ps.projector_radius()),
      nproj_(ps.num_projectors_lm()),
      matrix_({nproj_, sphere_.size()}),
      volume_element_(basis.volume_element()){

      std::vector<double> grid(sphere_.size()), proj(sphere_.size());

      // calculate the distance to the atom for each point
      for(int ipoint = 0; ipoint < sphere_.size(); ipoint++) grid[ipoint] = length(basis.rvector(sphere_.points()[ipoint]) - atom_position);
      
      for(int iproj = 0; iproj < ps.num_projectors_l(); iproj++){
				// interpolate the value of the radial part of the projectors to the sphere points
				ps.projector(iproj).value(sphere_.size(), grid, proj);
				
				int l = ps.projector_l(iproj);
				
				// now construct the projector with the spherical harmonics
				for(int m = -l; m <= l; m++){
					for(int ipoint = 0; ipoint < sphere_.size(); ipoint++){
						auto point = basis.rvector(sphere_.points()[ipoint]) - atom_position;
						matrix_[iproj][ipoint] = proj[ipoint]*math::spherical_harmonic(l, m, point);
					}
					
					kb_coeff_.push_back(ps.kb_coeff(iproj));
				}
				
      }
			
    }

    template <class array_dim4>
    void apply(const states::ks_states & st, const array_dim4 & phi, array_dim4 & vnlphi) const {

			//DATAOPERATIONS
			
      boost::multi::array<states::ks_states::coeff_type, 2> sphere_phi({sphere_.size(), st.num_states()});

      sphere_.gather(phi, sphere_phi);

      boost::multi::array<states::ks_states::coeff_type, 2> projections({nproj_, st.num_states()});

      //OPTIMIZATION: these two operations should be done by dgemm
      for(int iproj = 0; iproj < nproj_; iproj++){
				for(int ist = 0; ist < st.num_states(); ist++){
					complex aa = 0.0;
					for(int ipoint = 0; ipoint < sphere_.size(); ipoint++) aa += matrix_[iproj][ipoint]*sphere_phi[ipoint][ist];
					projections[iproj][ist] = aa*kb_coeff_[iproj]*volume_element_;
				}
      }

      for(int ipoint = 0; ipoint < sphere_.size(); ipoint++){
				for(int ist = 0; ist < st.num_states(); ist++){
					complex aa = 0.0;
					for(int iproj = 0; iproj < nproj_; iproj++) aa += matrix_[iproj][ipoint]*projections[iproj][ist];
					sphere_phi[ipoint][ist] = aa;
				}
      }
      
      sphere_.scatter_add(sphere_phi, vnlphi);
			
    }

    int num_projectors() const {
      return nproj_;
    }
    
  private:

    basis::spherical_grid sphere_;
    int nproj_;
    boost::multi::array<double, 2> matrix_;
    std::vector<double> kb_coeff_;
    double volume_element_;
    
  };
  
}

#endif

