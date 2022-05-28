/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SUBSPACE_DIAGONALIZATION
#define INQ__OPERATIONS__SUBSPACE_DIAGONALIZATION

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <inq_config.h>
#include <operations/diagonalize.hpp>
#include <operations/overlap.hpp>
#include <operations/rotate.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace ground_state {

template <class hamiltonian_type, class field_set_type>
auto subspace_diagonalization(const hamiltonian_type & ham, field_set_type & phi){
	CALI_CXX_MARK_FUNCTION;

	auto subspace_hamiltonian = operations::overlap(phi, ham(phi));
	auto eigenvalues = operations::diagonalize(subspace_hamiltonian);
	operations::rotate(subspace_hamiltonian, phi);
	return +eigenvalues({phi.set_part().start(), phi.set_part().end()});
}

}
}

#ifdef INQ_GROUND_STATE_SUBSPACE_DIAGONALIZATION_UNIT_TEST
#undef INQ_GROUND_STATE_SUBSPACE_DIAGONALIZATION_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::subspace_diagonalization", "[subspace_diagonalization]") {

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif
