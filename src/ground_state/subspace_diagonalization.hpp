/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SUBSPACE_DIAGONALIZATION
#define INQ__OPERATIONS__SUBSPACE_DIAGONALIZATION

/*
 Copyright (C) 2019 Xavier Andrade

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

#include <config.h>
#include <operations/overlap.hpp>
#include <operations/diagonalize.hpp>
#ifdef HAVE_CUDA
#include <multi/adaptors/blas/cuda.hpp> // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>

namespace inq {
namespace ground_state {

template <class hamiltonian_type, class field_set_type>
void subspace_diagonalization(const hamiltonian_type & ham, field_set_type & phi){
	
	auto subspace_hamiltonian = operations::overlap(phi, ham(phi));
	auto eigenvalues = operations::diagonalize(subspace_hamiltonian);
	
	//OPTIMIZATION: here we don't need to make a full copy. We can
	//divide into blocks over point index (second dimension of phi).
	using boost::multi::blas::gemm;
	using boost::multi::blas::transposed;
	
	phi.matrix() = gemm(phi.matrix(), transposed(subspace_hamiltonian));
	
}

}
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::subspace_diagonalization", "[subspace_diagonalization]") {

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif
