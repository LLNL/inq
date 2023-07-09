/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__CURRENT
#define INQ__OBSERVABLES__CURRENT

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Yifan Yao
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <operations/gradient.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <physics/constants.hpp>

namespace inq {
namespace observables {

template <typename HamiltonianType>
basis::field<basis::real_space, vector3<double, covariant>> current_density(const systems::ions & ions, systems::electrons const & electrons, HamiltonianType const & ham){

	basis::field<basis::real_space, vector3<double, covariant>> cdensity(electrons.density_basis());
	cdensity.fill(vector3<double, covariant>{0.0, 0.0, 0.0});

	auto iphi = 0;
	for(auto & phi : electrons.kpin()){
		
		auto gphi = operations::gradient(phi, /* factor = */ -1.0, /* shift = */ phi.kpoint() + ham.uniform_vector_potential());

		ham.projectors_all().position_commutator(phi, gphi, phi.kpoint() + ham.uniform_vector_potential());
		
    gpu::run(phi.basis().part().local_size(),
             [nst = phi.set_part().local_size(), occ = begin(electrons.occupations()[iphi]),
              ph = begin(phi.matrix()), gph = begin(gphi.matrix()), cdens = begin(cdensity.linear())] GPU_LAMBDA (auto ip){
               for(int ist = 0; ist < nst; ist++) cdens[ip] += 0.5*occ[ist]*imag(conj(ph[ip][ist])*gph[ip][ist] - conj(gph[ip][ist])*ph[ip][ist]);
             });
		iphi++;
	}
  
	cdensity.all_reduce(electrons.kpin_states_comm());
	return cdensity;
}

template <typename HamiltonianType>
auto current(const systems::ions & ions, systems::electrons const & electrons, HamiltonianType const & ham){
  return operations::integral(current_density(ions, electrons, ham));
}

}
}
#endif

#ifdef INQ_OBSERVABLES_CURRENT_UNIT_TEST
#undef INQ_OBSERVABLES_CURRENT_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;
		
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	auto par = input::parallelization(comm).states();

	{
		systems::ions ions(systems::cell::orthorhombic(6.0_b, 10.0_b, 6.0_b));
		systems::electrons electrons(par, ions, options::electrons{}.cutoff(15.0_Ha).extra_electrons(20.0));
		ground_state::initial_guess(ions, electrons);
		hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(), ions, 0.0, /* use_ace = */ true);
		
		SECTION("Gamma - no atoms"){

			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};
			
			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == -1.9018202326_a);
			CHECK(cur[1] == -5.8744437947_a);
			CHECK(cur[2] == -2.0854607655_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == -14.3800491961_a);
			CHECK(cur[1] ==  19.0820141323_a);
			CHECK(cur[2] == -39.5201476561_a);
			
		}

		SECTION("Gamma - no atoms - zero paramagnetic"){
			
			for(auto & phi : electrons.kpin()) phi.fill(1.0);
			
			auto charge = operations::integral(electrons.density());
			
			CHECK(charge == 20.0_a);
						
			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};

			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == Approx(0.0).margin(1e-12));
			CHECK(cur[1] == Approx(0.0).margin(1e-12));
			CHECK(cur[2] == Approx(0.0).margin(1e-12));
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);
			auto cur_an = -charge*ham.uniform_vector_potential()*ions.cell().volume();

			CHECK(cur[0] == Approx(cur_an[0]));
			CHECK(cur[1] == Approx(cur_an[1]));
			CHECK(cur[2] == Approx(cur_an[2]));						
			
		}
	}

	{
		systems::ions ions(systems::cell::orthorhombic(6.0_b, 10.0_b, 6.0_b));
		systems::electrons electrons(par, ions, options::electrons{}.cutoff(15.0_Ha).extra_electrons(20.0), input::kpoints::point(0.25, 0.25, 0.25));
		ground_state::initial_guess(ions, electrons);
		hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(), ions, 0.0, /* use_ace = */ true);
		
		SECTION("1/4 1/4 1/4 - no atoms"){
			
			auto cur = observables::current(ions, electrons, ham);

			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};
			
			CHECK(cur[0] == -21.5025764534_a);
			CHECK(cur[1] == -25.4752000156_a);
			CHECK(cur[2] == -21.6862169863_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == -33.9808054170_a);
			CHECK(cur[1] ==  -0.5187420885_a);
			CHECK(cur[2] == -59.1209038769_a);
			
		}

		SECTION("1/4 1/4 1/4 - no atoms - zero paramagnetic"){
			
			for(auto & phi : electrons.kpin()) phi.fill(1.0);
			
			auto charge = operations::integral(electrons.density());
			
			CHECK(charge == 20.0_a);

			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};
			
			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == Approx(-charge*0.25*2*M_PI*ions.cell().volume()));
			CHECK(cur[1] == Approx(-charge*0.25*2*M_PI*ions.cell().volume()));
			CHECK(cur[2] == Approx(-charge*0.25*2*M_PI*ions.cell().volume()));
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);
			auto cur_an = -charge*(vector3<double, covariant>{0.25, 0.25, 0.25}*2*M_PI + ham.uniform_vector_potential())*ions.cell().volume();

			CHECK(cur[0] == Approx(cur_an[0]));
			CHECK(cur[1] == Approx(cur_an[1]));
			CHECK(cur[2] == Approx(cur_an[2]));						
			
		}
	}

	
}
#endif
