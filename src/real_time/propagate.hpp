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
#include <observables/dipole.hpp>
#include <real_time/result.hpp>

namespace real_time {
	
	real_time::result propagate(systems::electrons & electrons, const input::interaction & inter){

		result res;
		
		const double dt = 0.055;
		
		const int numsteps = 100;
		
		auto density = density::calculate(electrons.states_.occupations(), electrons.phi_, electrons.density_basis_);
		
		hamiltonian::ks_hamiltonian<basis::real_space> ham(electrons.states_basis_, electrons.ions_.cell(), electrons.atomic_pot_, electrons.ions_.geo(), electrons.states_.num_states(), inter.exchange_coefficient());
		hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_);
		hamiltonian::energy energy;
		
		sc.update_ionic_fields(electrons.ions_, electrons.atomic_pot_);
		
		ham.scalar_potential = sc.ks_potential(density, energy);

		auto eigenvalues = operations::overlap_diagonal(electrons.phi_, ham(electrons.phi_));;
		energy.eigenvalues = operations::sum(electrons.states_.occupations(), eigenvalues, [](auto occ, auto ev){ return occ*real(ev); });
		
		tfm::format(std::cout, "step %9d :  t =  %9.3f e = %.12f\n", 0, 0.0, energy.total());

		res.time.push_back(0.0);
		res.energy.push_back(energy.total());
		res.dipole.push_back(observables::dipole(density));
		
		for(int istep = 1; istep <= numsteps; istep++){

			{
				//propagate half step and full step with H(t)
				auto fullstep_phi = operations::exponential_2_for_1(ham, complex(0.0, dt), complex(0.0, dt/2.0), electrons.phi_);
				
				//calculate H(t + dt) from the full step propagation
				density = density::calculate(electrons.states_.occupations(), fullstep_phi, electrons.density_basis_);
				ham.scalar_potential = sc.ks_potential(density, energy);
			}

			//propagate the other half step with H(t + dt)
			operations::exponential_in_place(ham, complex(0.0, dt/2.0), electrons.phi_);

			auto eigenvalues = operations::overlap_diagonal(electrons.phi_, ham(electrons.phi_));;
			energy.eigenvalues = operations::sum(electrons.states_.occupations(), eigenvalues, [](auto occ, auto ev){ return occ*real(ev); });
			
			tfm::format(std::cout, "step %9d :  t =  %9.3f e = %.12f\n", istep, istep*dt, energy.total());

			res.time.push_back(istep*dt);
			res.energy.push_back(energy.total());
			res.dipole.push_back(observables::dipole(density));			
			
		}

		return res;
	}
}

#endif
