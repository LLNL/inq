/* -*- indent-tabs-mode: t -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/energy.hpp>
#include <basis/field_set.hpp>
#include <operations/randomize.hpp>
#include <operations/integral.hpp>
#include <operations/orthogonalize.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <input/rt.hpp>
#include <input/scf.hpp>
#include <ions/interaction.hpp>
#include <real_time/result.hpp>

namespace systems {
	class electrons;
}
		
namespace ground_state {
	hamiltonian::energy calculate(systems::electrons & electrons, const input::interaction & inter = {}, const input::scf & solver = {});
}

namespace real_time {
	real_time::result propagate(systems::electrons & electrons, const input::interaction & inter = {}, const input::rt & options = {});
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
		friend real_time::result real_time::propagate(systems::electrons & electrons, const input::interaction & inter, const input::rt & options);
		
	private:
		
		const systems::ions & ions_;
		basis::real_space states_basis_;
		basis::real_space density_basis_;
		hamiltonian::atomic_potential atomic_pot_;
		states::ks_states states_;

	public: //temporary hack to be able to apply a kick from main
		
		basis::field_set<basis::real_space, complex> phi_;

	};  
  
}

#endif

