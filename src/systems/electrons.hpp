/* -*- indent-tabs-mode: t -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/self_consistency.hpp>
#include <hamiltonian/energy.hpp>
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
#include <input/config.hpp>
#include <functionals/lda.hpp>
#include <ions/interaction.hpp>

namespace systems {

  class electrons {

  public:

		enum class error { NO_ELECTRONS };
		
    electrons(const systems::ions & ions_arg, const input::basis arg_basis_input, const input::config & conf):
      ions_(ions_arg),
      rs_(ions_.cell(), arg_basis_input),
      atomic_pot_(ions_.geo().num_atoms(), ions_.geo().atoms()),
      states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons() + conf.excess_charge, conf.extra_states),
      ham_(rs_, ions_.cell(), atomic_pot_, ions_.geo()),
      phi_(rs_, states_.num_states()),
			sc_(conf.theory){

      rs_.info(std::cout);  
      states_.info(std::cout);

			if(atomic_pot_.num_electrons() + conf.excess_charge == 0) throw error::NO_ELECTRONS;
 
			ham_.info(std::cout);

      operations::randomize(phi_);
			operations::orthogonalization(phi_);
    }

    auto calculate_ground_state(const double ecutprec = 4.0) {

			//for the moment I am putting here some parameters that should be configurable. XA
			const double mixing = 0.1;
			//

			hamiltonian::energy energy;

			operations::preconditioner prec(ecutprec);
			
      double old_energy = DBL_MAX;

			auto vexternal = atomic_pot_.local_potential(rs_, ions_.cell(), ions_.geo());
			
			auto density = operations::calculate_density(states_.occupations(), phi_);

			ham_.scalar_potential = sc_.ks_potential(vexternal, density, atomic_pot_.ionic_density(rs_, ions_.cell(), ions_.geo()), energy);
			::ions::interaction_energy(atomic_pot_.range_separation(), ions_.cell(), ions_.geo(), energy.ion, energy.self);
																		 
      for(int ii = 0; ii < 1000; ii++){

				operations::subspace_diagonalization(ham_, phi_);
				
				solvers::steepest_descent(ham_, prec, phi_);

				density = operations::calculate_density(states_.occupations(), phi_);

				auto vks = sc_.ks_potential(vexternal, density, atomic_pot_.ionic_density(rs_, ions_.cell(), ions_.geo()), energy);
				
				auto eigenvalues = operations::overlap_diagonal(ham_(phi_), phi_);

				auto potdiff = operations::integral_absdiff(vks, ham_.scalar_potential)/abs(operations::integral(vks));
				
				//DATAOPERATIONS
				energy.eigenvalues = 0.0;
				for(int ii = 0; ii < states_.num_states(); ii++) energy.eigenvalues += states_.occupations()[ii]*real(eigenvalues[ii]);

				std::cout << "SCF iter " << ii << ":  e = " << std::scientific << energy.total()
									<< "  de = " << energy.total() - old_energy << "  dvks = " << potdiff << std::endl;

				for(int ii = 0; ii < states_.num_states(); ii++){
					std::cout << " state " << ii << ":  occ = " << states_.occupations()[ii] << "  evalue = " << real(eigenvalues[ii]) << std::endl;
				}
				std::cout << std::endl;
				
				if(fabs(energy.total() - old_energy) < 1e-5) break;
				
				old_energy = energy.total();

				//DATAOPERATIONS
				for(long ii = 0; ii < rs_.size(); ii++){
					ham_.scalar_potential[ii] = mixing*vks[ii] + (1.0 - mixing)*ham_.scalar_potential[ii];
				}
				
      }

			energy.print(std::cout);

			return energy;			
    }
    
  private:

		const systems::ions & ions_;
    basis::real_space rs_;
    hamiltonian::atomic_potential atomic_pot_;
    states::ks_states states_;
    hamiltonian::ks_hamiltonian<basis::real_space> ham_;      
    basis::field_set<basis::real_space, complex> phi_;
		hamiltonian::self_consistency sc_;

  };  
  
}

#endif

