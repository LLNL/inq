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
#include <operations/orthogonalize.hpp>
#include <operations/preconditioner.hpp>
#include <operations/integral.hpp>
#include <operations/subspace_diagonalization.hpp>
#include <density/calculate.hpp>
#include <density/normalize.hpp>
#include <mixers/linear.hpp>
#include <mixers/pulay.hpp>
#include <eigensolvers/conjugate_gradient.hpp>
#include <eigensolvers/steepest_descent.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <ions/interaction.hpp>
#include <input/scf.hpp>

namespace systems {

  class electrons {

  public:

		enum class error { NO_ELECTRONS };
		
    electrons(const systems::ions & ions_arg, const input::basis arg_basis_input, const input::config & conf):
      ions_(ions_arg),
      states_basis_(ions_.cell(), arg_basis_input),
			density_basis_(states_basis_.refine(arg_basis_input.density_factor())),
      atomic_pot_(ions_.geo().num_atoms(), ions_.geo().atoms()),
      states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons() + conf.excess_charge, conf.extra_states),
			phi_(states_basis_, states_.num_states()){
			
      states_basis_.info(std::cout);  
      states_.info(std::cout);

			if(atomic_pot_.num_electrons() + conf.excess_charge == 0) throw error::NO_ELECTRONS;
 

      operations::randomize(phi_);
			operations::orthogonalize(phi_);
    }

    auto calculate_ground_state(const input::interaction & inter, const input::scf & solver = {}){

			hamiltonian::ks_hamiltonian<basis::real_space> ham(states_basis_, ions_.cell(), atomic_pot_, ions_.geo(), states_.num_states(), inter.exchange_coefficient());

			ham.info(std::cout);

			hamiltonian::self_consistency sc(inter, states_basis_, density_basis_);

			hamiltonian::energy energy;

			operations::preconditioner prec;

			auto mixer = solvers::linear_mixer<double>(solver.mixing());
			//auto mixer = solvers::pulay_mixer<double>(2, solver.mixing(), states_basis.part().local_size());
			
      double old_energy = DBL_MAX;

			sc.update_ionic_fields(ions_, atomic_pot_);
	
			auto density = atomic_pot_.atomic_electronic_density(density_basis_, ions_.cell(), ions_.geo());
			density::normalize(density, states_.total_charge());
			std::cout << "Integral of the density = " << operations::integral(density) << std::endl;

			ham.scalar_potential = sc.ks_potential(density, energy);
			
			energy.ion = ::ions::interaction_energy(ions_.cell(), ions_.geo(), atomic_pot_);
			
			//DATAOPERATIONS STL FILL
			std::fill(ham.exchange.hf_occupations.begin(), ham.exchange.hf_occupations.end(), 0.0);

			ham.exchange.hf_orbitals = 0.0;
			
			int conv_count = 0;
      for(int iiter = 0; iiter < 1000; iiter++){

				operations::subspace_diagonalization(ham, phi_);

				{
					auto fphi = operations::space::to_fourier(std::move(phi_));

					switch(solver.eigensolver()){

					case input::scf::scf_eigensolver::STEEPEST_DESCENT:
						solvers::steepest_descent(ham, prec, fphi);
						break;
						
					case input::scf::scf_eigensolver::CONJUGATE_GRADIENT:
						eigensolver::conjugate_gradient(ham, prec, fphi);
						break;

					default:
						assert(false);
					}
					
					phi_ = operations::space::to_real(std::move(fphi));
				}

				//update the Hartree-Fock operator, mixing the new and old orbitals
				
				//DATAOPERATIONS LOOP 1D
				for(int ii = 0; ii < phi_.num_elements(); ii++){
					ham.exchange.hf_orbitals.data()[ii] = (1.0 - solver.mixing())*ham.exchange.hf_orbitals.data()[ii] + solver.mixing()*phi_.data()[ii];
				}
				
				//probably the occupations should be mixed too
				ham.exchange.hf_occupations = states_.occupations();

				if(inter.self_consistent() and solver.mix_density()) {
					auto new_density = density::calculate(states_.occupations(), phi_, density_basis_);
					mixer(density.linear(), new_density.linear());
					density::normalize(density, states_.total_charge());
				} else {
					density = density::calculate(states_.occupations(), phi_, density_basis_);
				}
				
				auto vks = sc.ks_potential(density, energy);

				if(inter.self_consistent() and solver.mix_potential()) {
					mixer(ham.scalar_potential.linear(), vks.linear());
				} else {
					ham.scalar_potential = std::move(vks);
				}
				
				// calculate the new energy and print
				{
					
					auto residual = ham(phi_);
					auto eigenvalues = operations::overlap_diagonal(phi_, residual);
					operations::shift(eigenvalues, phi_, residual, -1.0);
					
					auto normres = operations::overlap_diagonal(residual);
					auto nl_me = operations::overlap_diagonal(ham.non_local(phi_), phi_);
					auto exchange_me = operations::overlap_diagonal(ham.exchange(phi_), phi_);

					auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
					
					energy.eigenvalues = operations::sum(states_.occupations(), eigenvalues, energy_term);
					energy.nonlocal = operations::sum(states_.occupations(), nl_me, energy_term);
					energy.hf_exchange = operations::sum(states_.occupations(), exchange_me, energy_term);

					auto potdiff = operations::integral_absdiff(vks, ham.scalar_potential)/fabs(operations::integral(vks));
					
					tfm::format(std::cout, "SCF iter %d :  e = %.12f  de = %5.0e dvks = %5.0e\n",
											iiter, energy.total(), energy.eigenvalues - old_energy, potdiff);
					
					for(int istate = 0; istate < states_.num_states(); istate++){
						tfm::format(std::cout, " state %4d  occ = %4.3f  evalue = %18.12f  res = %5.0e\n",
												istate + 1, states_.occupations()[istate], real(eigenvalues[istate]), real(normres[istate]));
					}

				}
				
				if(fabs(energy.eigenvalues - old_energy) < solver.energy_tolerance()){
					conv_count++;
					if(conv_count > 2) break;
				} else {
					conv_count = 0;
				}
				
				old_energy = energy.eigenvalues;
				
      }

			energy.print(std::cout);

			return energy;			
    }
    
  private:

		const systems::ions & ions_;
    basis::real_space states_basis_;
		basis::real_space density_basis_;
    hamiltonian::atomic_potential atomic_pot_;
    states::ks_states states_;
		basis::field_set<basis::real_space, complex> phi_;


  };  
  
}

#endif

