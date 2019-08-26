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
		const int nst = phi.set_size();
		boost::multi::array<typename field_set_type::value_type, 1> tmp(nst);

		for(int ipoint = 0; ipoint < phi.basis().size(); ipoint++){
			for(int ii = 0; ii < nst; ii++){
				typename field_set_type::value_type aa = 0.0;
				for(int jj = 0; jj < nst; jj++) aa += conj(subspace_hamiltonian[ii][jj])*phi[ipoint][jj];
				tmp[ii] = aa;
			}
			for(int ii = 0; ii < nst; ii++) phi[ipoint][ii] = tmp[ii];
		}

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
