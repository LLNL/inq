/* -*- indent-tabs-mode: t -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

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
#include <ground_state/result.hpp>
#include <real_time/result.hpp>
#include <systems/ions.hpp>
#include <density/calculate.hpp>
#include <density/normalize.hpp>

#include <cfloat>

namespace inq {

namespace systems {
	class electrons;
}

namespace ground_state {
	ground_state::result calculate(const inq::systems::ions & ions, inq::systems::electrons & electrons, const input::interaction & inter = {}, const input::scf & solver = {});
}

namespace real_time {
	real_time::result propagate(inq::systems::ions & ions, inq::systems::electrons & electrons, const input::interaction & inter = {}, const input::rt & options = {});
}

namespace systems {

	class electrons {
		
	public:
	
		enum class error { NO_ELECTRONS };

		electrons(boost::mpi3::communicator & comm, const inq::systems::ions & ions, const input::basis arg_basis_input, const input::config & conf = {}):
			full_comm_(comm, {}),
			states_comm_(full_comm_.axis(0)),
			basis_comm_(full_comm_.axis(1)),
			states_basis_(ions.cell(), arg_basis_input),
			density_basis_(states_basis_.refine(arg_basis_input.density_factor())),
			atomic_pot_(ions.geo().num_atoms(), ions.geo().atoms(), states_basis_.gcutoff()),
			states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons() + conf.excess_charge, conf.extra_states),
			phi_(states_basis_, states_.num_states()),
			density_(atomic_pot_.atomic_electronic_density(density_basis_, ions.cell(), ions.geo()))
		{

			if(comm.root()){
				states_basis_.info(std::cout);  
				states_.info(std::cout);
			}
			
			if(atomic_pot_.num_electrons() + conf.excess_charge == 0) throw error::NO_ELECTRONS;
 

			operations::randomize(phi_);
			operations::orthogonalize(phi_);
			
			density::normalize(density_, states_.total_charge());

			if(comm.root()) std::cout << "Integral of the density = " << operations::integral(density_) << std::endl;
			
    }


	public: //temporary hack to be able to apply a kick from main and avoid a bug in nvcc

		mutable boost::mpi3::cartesian_communicator<2> full_comm_;
		mutable boost::mpi3::cartesian_communicator<1> states_comm_;
		mutable boost::mpi3::cartesian_communicator<1> basis_comm_;
		basis::real_space states_basis_;
		basis::real_space density_basis_;
		hamiltonian::atomic_potential atomic_pot_;
		states::ks_states states_;
		basis::field_set<basis::real_space, complex> phi_;
		basis::field<basis::real_space, double> density_;

	};
  
}
}

#endif

