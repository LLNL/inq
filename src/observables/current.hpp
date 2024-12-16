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
	auto par = input::parallelization(comm);

	{
		systems::ions ions(systems::cell::orthorhombic(6.0_b, 10.0_b, 6.0_b));
		systems::electrons electrons(par, ions, options::electrons{}.cutoff(15.0_Ha).extra_electrons(20.0));
		ground_state::initial_guess(ions, electrons);
		hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(), ions, 0.0, /* use_ace = */ true);
		
		SECTION("Gamma - no atoms"){

			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};
			
			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == 62.2452955214_a);
			CHECK(cur[1] == -1.0723045428_a);
			CHECK(cur[2] == 55.1882949624_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] ==  42.2452955214_a);
			CHECK(cur[1] ==  38.9276954572_a);
			CHECK(cur[2] ==  -4.8117050376_a);
			
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
			
			CHECK(cur[0] ==  30.8293689855_a);
			CHECK(cur[1] == -32.4882310787_a);
			CHECK(cur[2] ==  23.7723684265_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] ==  10.8293689855_a);
			CHECK(cur[1] ==   7.5117689213_a);
			CHECK(cur[2] == -36.2276315735_a);
			
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

	{
		systems::ions ions(systems::cell::orthorhombic(6.0_b, 10.0_b, 6.0_b));
		ions.insert("Cu", {0.0_b, 0.0_b, 0.0_b});
		ions.insert("Ag", {2.0_b, 0.7_b, 0.0_b});	
		systems::electrons electrons(par, ions, options::electrons{}.cutoff(15.0_Ha));
		ground_state::initial_guess(ions, electrons);
		hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(), ions, 0.0, /* use_ace = */ true);
		
		SECTION("Gamma - atoms"){

			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};
			
			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == 125.0726107094_a);
			CHECK(cur[1] ==  4.7292765308_a);
			CHECK(cur[2] == 116.5420097097_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] ==  87.0638041144_a);
			CHECK(cur[1] ==  80.6922069845_a);
			CHECK(cur[2] ==   2.5258837095_a);
			
		}

		SECTION("Gamma - atoms - zero paramagnetic"){
			
			for(auto & phi : electrons.kpin()) phi.fill(1.0);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};

			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == Approx(0.0).margin(1e-12));
			CHECK(cur[1] == Approx(0.0).margin(1e-12));
			CHECK(cur[2] == Approx(0.0).margin(1e-12));
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == -13655.6186604198_a);
			CHECK(cur[1] ==  27311.2419971710_a);
			CHECK(cur[2] == -40966.6831884187_a);						
			
		}
	}

	{
		systems::ions ions(systems::cell::orthorhombic(6.0_b, 10.0_b, 6.0_b));
		ions.insert("Cu", {0.0_b, 0.0_b, 0.0_b});
		ions.insert("Ag", {2.0_b, 0.7_b, 0.0_b});	
		systems::electrons electrons(par, ions, options::electrons{}.cutoff(15.0_Ha), input::kpoints::point(0.25, 0.25, 0.25));
		ground_state::initial_guess(ions, electrons);
		hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(), ions, 0.0, /* use_ace = */ true);
		
		SECTION("Gamma - atoms"){

			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};
			
			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] ==  65.3579017395_a);
			CHECK(cur[1] == -54.9266021369_a);
			CHECK(cur[2] ==  56.8453806312_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] ==  27.3599465174_a);
			CHECK(cur[1] ==  21.0517395749_a);
			CHECK(cur[2] == -57.1448517995_a);
			
		}

		SECTION("Gamma - atoms - zero paramagnetic"){
			
			for(auto & phi : electrons.kpin()) phi.fill(1.0);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{0.0, 0.0, 0.0};

			auto cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == -21452.9748736864_a);
			CHECK(cur[1] == -21453.0100725570_a);
			CHECK(cur[2] == -21452.8830570262_a);
			
			ham.uniform_vector_potential() = vector3<double, covariant>{1.0, -2.0, 3.0};

			cur = observables::current(ions, electrons, ham);

			CHECK(cur[0] == -35093.1574770241_a);
			CHECK(cur[1] ==   5859.0105019392_a);
			CHECK(cur[2] == -62394.5776586859_a);
			
		}
	}

	
}
#endif
