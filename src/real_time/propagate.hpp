/* -*- indent-tabs-mode: t -*- */

#ifndef REAL_TIME__PROPAGATE
#define REAL_TIME__PROPAGATE

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
#include <operations/exponential.hpp>
#include <density/calculate.hpp>
#include <density/normalize.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <ions/interaction.hpp>
#include <input/scf.hpp>
#include <systems/electrons.hpp>

namespace real_time {
	
	void propagate(systems::electrons & electrons, const input::interaction & inter){
		
		const double dt = 0.01;
		
		const int numsteps = 100;
		
		auto density = density::calculate(electrons.states_.occupations(), electrons.phi_, electrons.density_basis_);
		
		hamiltonian::ks_hamiltonian<basis::real_space> ham(electrons.states_basis_, electrons.ions_.cell(), electrons.atomic_pot_, electrons.ions_.geo(), electrons.states_.num_states(), inter.exchange_coefficient());
		hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_);
		hamiltonian::energy energy;
		
		sc.update_ionic_fields(electrons.ions_, electrons.atomic_pot_);
		
		ham.scalar_potential = sc.ks_potential(density, energy);

		for(int istep = 0; istep < numsteps; istep++){
			
			electrons.phi_ = operations::exponential(ham, complex(0.0, dt), electrons.phi_);	

			ham.scalar_potential = sc.ks_potential(density, energy);

			auto eigenvalues = operations::overlap_diagonal(electrons.phi_, ham(electrons.phi_));;
			energy.eigenvalues = operations::sum(electrons.states_.occupations(), eigenvalues, [](auto occ, auto ev){ return occ*real(ev); });
			
			tfm::format(std::cout, "step %d :  e = %.12f\n", istep, energy.total());
			
		}
	}
}

#endif

