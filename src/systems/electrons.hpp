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
	class electrons;
}

namespace ground_state{
	hamiltonian::energy calculate(systems::electrons & electrons, const input::interaction & inter, const input::scf & solver = {});
}

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

		friend hamiltonian::energy ground_state::calculate(systems::electrons & electrons, const input::interaction & inter, const input::scf & solver);

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

