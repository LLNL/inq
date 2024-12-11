/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__INITIAL_GUESS
#define INQ__GROUND_STATE__INITIAL_GUESS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cfloat>

#include <systems/ions.hpp>
#include <operations/randomize.hpp>
#include <operations/orthogonalize.hpp>
#include <observables/density.hpp>
#include <systems/electrons.hpp>

namespace inq {
namespace ground_state {
	
void initial_guess(const systems::ions & ions, systems::electrons & electrons, std::optional<vector3<double>> const & magnet_dir = {}){

	int iphi = 0;
	for(auto & phi : electrons.kpin()) {
		operations::randomize(phi, iphi + electrons.kpin_part().start());
		operations::orthogonalize(phi);
		for(long ist = 0; ist < phi.local_spinor_set_size(); ist++) electrons.eigenvalues()[iphi][ist] = ist + phi.spinor_set_part().start() + (iphi + electrons.kpin_part().start())/double(electrons.kpin_part().size());

		iphi++;
	}
	
	electrons.update_occupations(electrons.eigenvalues());
	
	if(ions.size() > 0){
		electrons.spin_density() = electrons.atomic_pot().atomic_electronic_density(electrons.states_comm(), electrons.density_basis(), ions, electrons.states());
	} else {
		electrons.spin_density() = observables::density::calculate(electrons);
	}

	assert(fabs(operations::integral_sum(electrons.spin_density())) > 1e-16);
	
  observables::density::normalize(electrons.spin_density(), electrons.states().num_electrons());
  if (magnet_dir) {
	assert(electrons.spin_density().set_size() > 1);
	auto magnet_dir_ = {magnet_dir.value()[0], magnet_dir.value()[1], magnet_dir.value()[2]};
	observables::density::rotate_total_magnetization(electrons.spin_density(), magnet_dir_);
  }

}
}
}
#endif

#ifdef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST
#undef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	SECTION("Spin unpolarized initialization") {
		auto par = input::parallelization(comm);
		auto ions = systems::ions(systems::cell::cubic(10.0_b));
		ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
		auto electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_unpolarized());
		ground_state::initial_guess(ions, electrons);
		CHECK(Approx(operations::integral_sum(electrons.spin_density())).epsilon(1.e-10)     == 1.0);
	}

	SECTION("Spin polarized initialization") {
		auto par = input::parallelization(comm);
		auto ions = systems::ions(systems::cell::cubic(10.0_b));
		ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
		auto electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
		ground_state::initial_guess(ions, electrons);
		CHECK(Approx(operations::integral_sum(electrons.spin_density())).epsilon(1.e-10)     == 1.0);

		vector3 mag_dir = {0.0, 0.0, 1.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
		ground_state::initial_guess(ions, electrons, mag_dir);
		auto mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral_sum(electrons.spin_density())).epsilon(1.e-10)     == 1.0);
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                 == 1.0);
		CHECK(Approx(sqrt(mag[0]*mag[0]+mag[1]*mag[1])/sqrt(norm(mag))).epsilon(1.e-10)      == 0.0);

		mag_dir = {0.0, 0.0, -1.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
		ground_state::initial_guess(ions, electrons, mag_dir);
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral_sum(electrons.spin_density())).epsilon(1.e-10)     == 1.0);
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                 == -1.0);
		CHECK(Approx(sqrt(mag[0]*mag[0]+mag[1]*mag[1])/sqrt(norm(mag))).epsilon(1.e-10)      == 0.0);
	}

	SECTION("Spin non collinear initialization") {
		auto par = input::parallelization(comm);
		auto ions = systems::ions(systems::cell::cubic(10.0_b));
		ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
		auto electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons);
		CHECK(Approx(operations::integral_sum(electrons.spin_density())).epsilon(1.e-10)     == 1.0);

		vector3 mag_dir = {0.0, 0.0, 1.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		auto ch_density = observables::density::total(electrons.spin_density());
		auto mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                       == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                 == 0.0);
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                 == 0.0);
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                 == 1.0);

		mag_dir = {0.0, 0.0, -1.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		ch_density = observables::density::total(electrons.spin_density());
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                       == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                 == 0.0);
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                 == 0.0);
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                 == -1.0);

		mag_dir = {1.0, 1.0, 0.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		ch_density = observables::density::total(electrons.spin_density());
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                        == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(2.0));
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(2.0));
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 0.0);

		mag_dir = {-1.0, 1.0, 0.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		ch_density = observables::density::total(electrons.spin_density());
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                        == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                  == -1.0/sqrt(2.0));
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(2.0));
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 0.0);

		mag_dir = {1.0, -1.0, 0.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		ch_density = observables::density::total(electrons.spin_density());
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                        == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(2.0));
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                  == -1.0/sqrt(2.0));
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 0.0);

		mag_dir = {-1.0, -1.0, 0.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		ch_density = observables::density::total(electrons.spin_density());
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                        == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                  == -1.0/sqrt(2.0));
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                  == -1.0/sqrt(2.0));
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 0.0);

		mag_dir = {1.0, 1.0, 1.0};
		electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
		ground_state::initial_guess(ions, electrons, mag_dir);
		ch_density = observables::density::total(electrons.spin_density());
		mag = observables::total_magnetization(electrons.spin_density());
		CHECK(Approx(operations::integral(ch_density)).epsilon(1.e-10)                        == 1.0);
		CHECK(Approx(mag[0]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(3.0));
		CHECK(Approx(mag[1]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(3.0));
		CHECK(Approx(mag[2]/sqrt(norm(mag))).epsilon(1.e-10)                                  == 1.0/sqrt(3.0));
	}
}
#endif

