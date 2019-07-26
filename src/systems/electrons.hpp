/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <basis/field_set.hpp>
#include <operations/randomize.hpp>
#include <operations/overlap.hpp>
#include <operations/scal.hpp>
#include <operations/orthogonalization.hpp>
#include <operations/preconditioner.hpp>
#include <operations/calculate_density.hpp>
#include <operations/integral.hpp>
#include <solvers/steepest_descent.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <functionals/lda.hpp>

namespace systems {

  class electrons {

  public:
    
    electrons(const systems::ions & ions_arg, const input::basis arg_basis_input):
      ions_(ions_arg),
      rs_(ions_.cell(), arg_basis_input),
      atomic_pot_(ions_.geo().num_atoms(), ions_.geo().atoms()),
      states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons()),
      ham_(rs_, ions_.cell(), atomic_pot_, ions_.geo()),
      phi_(rs_, states_.num_states()),
			prec_(4){

      rs_.info(std::cout);  
      states_.info(std::cout);
      ham_.info(std::cout);

      operations::randomize(phi_);
    }

    void calculate_ground_state() {

			solvers::poisson<basis::real_space> poisson_solver;
			
      double old_energy = DBL_MAX;

			auto density = operations::calculate_density(states_.occupations(), phi_);
			
      for(int ii = 0; ii < 1000; ii++){

				auto hphi = ham_(phi_);
				
				auto overlap = operations::overlap_diagonal(hphi, phi_);
				
				//DATAOPERATIONS
				double energy = 0.0;
				for(int ii = 0; ii < states_.num_states(); ii++) energy += real(overlap[ii]);
				
				std::cout << ii << '\t' << std::scientific << energy << '\t' << energy - old_energy << std::endl;
				
				if(fabs(energy - old_energy) < 1e-7) break;
				
				old_energy = energy;
				
				solvers::steepest_descent(states_, ham_, prec_, phi_);

				density = operations::calculate_density(states_.occupations(), phi_);

				auto vhartree = poisson_solver(density);

				basis::field<basis::real_space, double> exc(vhartree.basis());
				basis::field<basis::real_space, double> vxc(vhartree.basis());

				functionals::lda::xc_unpolarized(density.basis().size(), density, exc, vxc);
				
      }
    }

    auto calculate_energy() {
      
      operations::scal_invsqrt(operations::overlap_diagonal(phi_), phi_);
      
      auto hphi = ham_(phi_);
      
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
    basis::field_set<basis::real_space, complex> phi_;
		operations::preconditioner prec_;

  };  
  
}

#endif

