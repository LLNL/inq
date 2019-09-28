/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SUBSPACE_DIAGONALIZATION
#define OPERATIONS__SUBSPACE_DIAGONALIZATION

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


namespace operations {

	template <class hamiltonian_type, class field_set_type>
  void subspace_diagonalization(const hamiltonian_type & ham, field_set_type & phi){

		auto subspace_hamiltonian = overlap(phi, ham(phi));
		auto eigenvalues = diagonalize(subspace_hamiltonian);

		//this should be done by multi + blas
		//DATAOPERATIONS

		//OPTIMIZATION: here we don't need to make a full copy. We can
		//divide into blocks over point index (second dimension of phi).
		boost::multi::array<typename field_set_type::value_type, 2> tmp = phi;
		
		boost::multi::blas::gemm('T', 'N', 1.0, subspace_hamiltonian, tmp, 0.0, phi);
 
  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::subspace_diagonalization", "[subspace_diagonalization]") {

	using namespace Catch::literals;
	
}

#endif

#endif
