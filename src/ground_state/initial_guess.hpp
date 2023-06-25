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
	
void initial_guess(const systems::ions & ions, systems::electrons & electrons){

	int iphi = 0;
	for(auto & phi : electrons.kpin()) {
		operations::randomize(phi, iphi + electrons.kpin_part().start());
		operations::orthogonalize(phi);
		for(long ist = 0; ist < phi.local_set_size(); ist++) electrons.eigenvalues()[iphi][ist] = ist + phi.set_part().start() + (iphi + electrons.kpin_part().start())/double(electrons.kpin_size());

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

}
}
}
#endif

#ifdef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST
#undef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif

