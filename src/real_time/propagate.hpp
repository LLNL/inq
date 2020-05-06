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
		
		hamiltonian::ks_hamiltonian<basis::real_space> ham(electrons.states_basis_, electrons.ions_.cell(), electrons.atomic_pot_, electrons.ions_.geo(), electrons.states_.num_states(), inter.exchange_coefficient());
		
		hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_);
		
		sc.update_ionic_fields(electrons.ions_, electrons.atomic_pot_);

		auto density = density::calculate(electrons.states_.occupations(), electrons.phi_, electrons.density_basis_);

	}
	
}

#endif

