/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SUBSPACE_DIAGONALIZATION
#define INQ__OPERATIONS__SUBSPACE_DIAGONALIZATION

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <matrix/diagonalize.hpp>
#include <operations/overlap.hpp>
#include <operations/rotate.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace ground_state {

template <class hamiltonian_type, class field_set_type>
auto subspace_diagonalization(const hamiltonian_type & ham, field_set_type & phi){
	CALI_CXX_MARK_FUNCTION;

	auto subspace_hamiltonian = operations::overlap(phi, ham(phi));
	auto eigenvalues = matrix::diagonalize(subspace_hamiltonian);
	operations::rotate(subspace_hamiltonian, phi);
	return +eigenvalues({phi.spinor_set_part().start(), phi.spinor_set_part().end()});
}

}
}
#endif

#ifdef INQ_GROUND_STATE_SUBSPACE_DIAGONALIZATION_UNIT_TEST
#undef INQ_GROUND_STATE_SUBSPACE_DIAGONALIZATION_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <operations/randomize.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif
