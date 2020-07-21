/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__CALCULATE
#define INQ__GROUND_STATE__CALCULATE

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
#include <density/calculate.hpp>
#include <density/normalize.hpp>
#include <mixers/linear.hpp>
#include <mixers/pulay.hpp>
#include <mixers/broyden.hpp>
#include <eigensolvers/conjugate_gradient.hpp>
#include <eigensolvers/steepest_descent.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <ions/interaction.hpp>
#include <input/scf.hpp>
#include <observables/dipole.hpp>
#include <systems/electrons.hpp>
#include <ground_state/result.hpp>
#include <ground_state/subspace_diagonalization.hpp>

namespace inq {
namespace ground_state {
	
	ground_state::result calculate(const systems::ions & ions, systems::electrons & electrons, const input::interaction & inter, const input::scf & solver){
		
		hamiltonian::ks_hamiltonian<basis::real_space> ham(electrons.states_basis_, ions.cell(), electrons.atomic_pot_, inter.fourier_pseudo_value(), ions.geo(), electrons.states_.num_states(), inter.exchange_coefficient());
		
		ham.info(std::cout);
		
		hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_);

		ground_state::result res;
		
		operations::preconditioner prec;

		mixers::base<double> * mixer = nullptr;

		switch(solver.mixing_algorithm()){
		case input::scf::mixing_algo::LINEAR:
			mixer = new mixers::linear<double>(solver.mixing());
			break;
		case input::scf::mixing_algo::PULAY:
			mixer = new mixers::pulay<double>(4, solver.mixing(), electrons.states_basis_.part().local_size());
			break;
		case input::scf::mixing_algo::BROYDEN:
			mixer = new mixers::broyden<double>(4, solver.mixing(), electrons.states_basis_.part().local_size());
			break;
		}
		
		double old_energy = DBL_MAX;
		
		sc.update_ionic_fields(ions, electrons.atomic_pot_);
		
		ham.scalar_potential = sc.ks_potential(electrons.density_, res.energy);
		
		res.energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);
		
		//DATAOPERATIONS STL FILL
		std::fill(ham.exchange.hf_occupations.begin(), ham.exchange.hf_occupations.end(), 0.0);
		
		ham.exchange.hf_orbitals = 0.0;
		
		int conv_count = 0;
		for(int iiter = 0; iiter < 1000; iiter++){
			
			subspace_diagonalization(ham, electrons.phi_);
			
			{
				auto fphi = operations::space::to_fourier(std::move(electrons.phi_));
				
				switch(solver.eigensolver()){
					
				case input::scf::scf_eigensolver::STEEPEST_DESCENT:
					eigensolvers::steepest_descent(ham, prec, fphi);
					break;
					
				case input::scf::scf_eigensolver::CONJUGATE_GRADIENT:
					eigensolvers::conjugate_gradient(ham, prec, fphi);
					break;
					
				default:
					assert(false);
				}
				
				electrons.phi_ = operations::space::to_real(std::move(fphi));
				
			}
			
			//update the Hartree-Fock operator, mixing the new and old orbitals
			
			//DATAOPERATIONS LOOP 1D
			for(int ii = 0; ii < electrons.phi_.num_elements(); ii++){
				ham.exchange.hf_orbitals.data()[ii] = (1.0 - solver.mixing())*ham.exchange.hf_orbitals.data()[ii] + solver.mixing()*electrons.phi_.data()[ii];
			}
			
			//probably the occupations should be mixed too
			ham.exchange.hf_occupations = electrons.states_.occupations();
			
			if(inter.self_consistent() and solver.mix_density()) {
				auto new_density = density::calculate(electrons.states_.occupations(), electrons.phi_, electrons.density_basis_);
				mixer->operator()(electrons.density_.linear(), new_density.linear());
				density::normalize(electrons.density_, electrons.states_.total_charge());
			} else {
				electrons.density_ = density::calculate(electrons.states_.occupations(), electrons.phi_, electrons.density_basis_);
			}
			
			auto vks = sc.ks_potential(electrons.density_, res.energy);
			
			if(inter.self_consistent() and solver.mix_potential()) {
				mixer->operator()(ham.scalar_potential.linear(), vks.linear());
			} else {
				ham.scalar_potential = std::move(vks);
			}
			
			// calculate the new energy and print
			{
				
				auto residual = ham(electrons.phi_);
				auto eigenvalues = operations::overlap_diagonal(electrons.phi_, residual);
				operations::shift(-1.0, eigenvalues, electrons.phi_, residual);
				
				auto normres = operations::overlap_diagonal(residual);
				auto nl_me = operations::overlap_diagonal(ham.non_local(electrons.phi_), electrons.phi_);
				auto exchange_me = operations::overlap_diagonal(ham.exchange(electrons.phi_), electrons.phi_);
				
				auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
				
				res.energy.eigenvalues = operations::sum(electrons.states_.occupations(), eigenvalues, energy_term);
				res.energy.nonlocal = operations::sum(electrons.states_.occupations(), nl_me, energy_term);
				res.energy.hf_exchange = operations::sum(electrons.states_.occupations(), exchange_me, energy_term);
				
				auto potdiff = operations::integral_absdiff(vks, ham.scalar_potential)/fabs(operations::integral(vks));

				if(solver.verbose_output()){
					
					tfm::format(std::cout, "SCF iter %d :  e = %.12f  de = %5.0e dvks = %5.0e\n",
											iiter, res.energy.total(), res.energy.eigenvalues - old_energy, potdiff);
					
					for(int istate = 0; istate < electrons.states_.num_states(); istate++){
						tfm::format(std::cout, " state %4d  occ = %4.3f  evalue = %18.12f  res = %5.0e\n",
												istate + 1, electrons.states_.occupations()[istate], real(eigenvalues[istate]), real(normres[istate]));
					}
				}
				
			}
			
			if(fabs(res.energy.eigenvalues - old_energy) < solver.energy_tolerance()){
				conv_count++;
				if(conv_count > 2) break;
			} else {
				conv_count = 0;
			}
			
			old_energy = res.energy.eigenvalues;
			
		}

		delete mixer;

		if(solver.verbose_output()){
			res.energy.print(std::cout);
		}

		if(ions.cell().periodic_dimensions() == 0){
			res.dipole = observables::dipole(ions, electrons);
		} else {
			res.dipole = 0.0;
		}

		//make sure we have a density consistet with phi
		electrons.density_ = density::calculate(electrons.states_.occupations(), electrons.phi_, electrons.density_basis_);

		return res;
	}
}
}

#endif

