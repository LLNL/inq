/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <basis/coefficients_set.hpp>
#include <operations/randomize.hpp>
#include <operations/overlap.hpp>
#include <operations/scal.hpp>
#include <operations/orthogonalization.hpp>
#include <solvers/steepest_descent.hpp>
#include <math/complex.hpp>

namespace systems {

  class electrons {

  public:
    
    electrons(const systems::ions & ions_arg, const double & ecut):
      ions_(ions_arg),
      rs_(ions_.cell(), ecut),
      atomic_pot_(ions_.geo().num_atoms(), ions_.geo().atoms()),
      states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons()),
      ham_(rs_, ions_.cell(), atomic_pot_, ions_.geo()),
      phi_(rs_, states_.num_states()){

      rs_.info(std::cout);  
      states_.info(std::cout);
      ham_.info(std::cout);

      operations::randomize(phi_);
    }

    void calculate_ground_state() {

      double old_energy = DBL_MAX;
      
      for(int ii = 0; ii < 2000; ii++){
				
				operations::scal_invsqrt(operations::overlap_diagonal(phi_), phi_);
				
				auto hphi = ham_(states_, phi_);
				
				auto overlap = operations::overlap_diagonal(hphi, phi_);
				
				//DATAOPERATIONS
				double energy = 0.0;
				for(int ii = 0; ii < states_.num_states(); ii++) energy += real(overlap[ii]);
				
				std::cout << ii << '\t' << std::scientific << energy << std::endl;
				
				if(fabs(energy - old_energy) < 1e-7) break;
				
				old_energy = energy;
				
				solvers::steepest_descent(states_, ham_, phi_);
				
      }
    }

    auto calculate_energy() {
      
      operations::scal_invsqrt(operations::overlap_diagonal(phi_), phi_);
      
      auto hphi = ham_(states_, phi_);
      
      auto overlap = operations::overlap_diagonal(hphi, phi_);

      //DATAOPERATIONS
      double energy = 0.0;
      for(int ii = 0; ii < states_.num_states(); ii++) energy += real(overlap[ii]);
      return energy;
    }


    
  private:

		const systems::ions & ions_;
    basis::real_space rs_;
    hamiltonian::atomic_potential atomic_pot_;
    states::ks_states states_;
    hamiltonian::ks_hamiltonian<basis::real_space> ham_;      
    basis::coefficients_set<basis::real_space, complex> phi_;
  

  };  
  
}

#endif

