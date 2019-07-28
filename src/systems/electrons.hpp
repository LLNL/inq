/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/ks_potential.hpp>
#include <basis/field_set.hpp>
#include <operations/randomize.hpp>
#include <operations/overlap.hpp>
#include <operations/scal.hpp>
#include <operations/orthogonalization.hpp>
#include <operations/preconditioner.hpp>
#include <operations/calculate_density.hpp>
#include <operations/integral.hpp>
#include <operations/subspace_diagonalization.hpp>
#include <solvers/steepest_descent.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <functionals/lda.hpp>

namespace systems {

  class electrons {

  public:

		enum class error { NO_ELECTRONS };
		
    electrons(const systems::ions & ions_arg, const input::basis arg_basis_input, const int extra_states = 0, const double excess_charge = 0.0):
      ions_(ions_arg),
      rs_(ions_.cell(), arg_basis_input),
      atomic_pot_(ions_.geo().num_atoms(), ions_.geo().atoms()),
      states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons() + excess_charge, extra_states),
      ham_(rs_, ions_.cell(), atomic_pot_, ions_.geo()),
      phi_(rs_, states_.num_states()){

      rs_.info(std::cout);  
      states_.info(std::cout);

			if(atomic_pot_.num_electrons()  + excess_charge == 0){
				throw error::NO_ELECTRONS;
			}

			ham_.info(std::cout);

      operations::randomize(phi_);
    }

    void calculate_ground_state() {

			//for the moment I am putting here some parameters that should be configurable. XA
			const double ecutprec = 4.0;			
			const double mixing = 0.1;
			//

			double energy, eeigenvalues, eexternal, ehartree, exc, intvxc;
			
			operations::preconditioner prec(ecutprec);
			
      double old_energy = DBL_MAX;

			auto vexternal = atomic_pot_.local_potential(rs_, ions_.cell(), ions_.geo());
			auto density = operations::calculate_density(states_.occupations(), phi_);

			ham_.scalar_potential = hamiltonian::ks_potential(vexternal, density, eexternal, ehartree, exc, intvxc);
			
      for(int ii = 0; ii < 1000; ii++){

				operations::subspace_diagonalization(ham_, phi_);
				
				solvers::steepest_descent(states_, ham_, prec, phi_);

				density = operations::calculate_density(states_.occupations(), phi_);

				auto vks = hamiltonian::ks_potential(vexternal, density, eexternal, ehartree, exc, intvxc);
				
				auto eigenvalues = operations::overlap_diagonal(ham_(phi_), phi_);

				auto potdiff = operations::integral_absdiff(vks, ham_.scalar_potential)/abs(operations::integral(vks));
				
				//DATAOPERATIONS
				eeigenvalues = 0.0;
				for(int ii = 0; ii < states_.num_states(); ii++) eeigenvalues += states_.occupations()[ii]*real(eigenvalues[ii]);
				energy = eeigenvalues + eexternal + ehartree + exc - intvxc;
				
				std::cout << "SCF iter " << ii << ":  e = " << std::scientific << energy << "  de = " << energy - old_energy << "  dvks = " << potdiff << std::endl;

				for(int ii = 0; ii < states_.num_states(); ii++){
					std::cout << " state " << ii << ":  occ = " << states_.occupations()[ii] << "  evalue = " << real(eigenvalues[ii]) << std::endl;
				}
				std::cout << std::endl;
				
				if(fabs(energy - old_energy) < 1e-5) break;
				
				old_energy = energy;

				//DATAOPERATIONS
				for(long ii = 0; ii < rs_.size(); ii++){
					ham_.scalar_potential[ii] = mixing*vks[ii] + (1.0 - mixing)*ham_.scalar_potential[ii];
				}
				
      }

			std::cout << std::endl;
			std::cout << "  total       = " << energy       << std::endl;
			std::cout << "  eigenvalues = " << eeigenvalues << std::endl;
			std::cout << "  external    = " << eexternal    << std::endl;
			std::cout << "  hartree     = " << ehartree     << std::endl;
			std::cout << "  xc          = " << exc          << std::endl;
			std::cout << "  intnvxc     = " << intvxc       << std::endl;
			std::cout << std::endl;
			
    }
    
  private:

		const systems::ions & ions_;
    basis::real_space rs_;
    hamiltonian::atomic_potential atomic_pot_;
    states::ks_states states_;
    hamiltonian::ks_hamiltonian<basis::real_space> ham_;      
    basis::field_set<basis::real_space, complex> phi_;

  };  
  
}

#endif

