/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SOLVERS_STEEPEST_DESCENT
#define SOLVERS_STEEPEST_DESCENT

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

#include <math/complex.hpp>
#include <math/d3vector.hpp>
#include <multi/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/shift.hpp>

namespace solvers {

	template <class operator_type>
	void steepest_descent(const states::ks_states st, const basis::plane_wave & basis, const operator_type & ham, states::coefficients & phi){

		//calculate the residual
		
		auto residual = ham.apply(st, basis, phi);

		auto eigenvalues = operations::overlap_diagonal(st, basis, residual, phi);
		auto norm = operations::overlap_diagonal(st, basis, phi);

		auto lambda(eigenvalues);
		
		//DATAOPERATIONS
		for(int ist = 0; ist < st.num_states(); ist++) lambda[ist] /= -norm[ist];

		operations::shift(st, basis, lambda, phi, residual);

		//OPTIMIZATIONS: precondition the residual here


		//now calculate the step size
		auto hresidual = ham.apply(st, basis, residual);

		auto m1 = operations::overlap_diagonal(st, basis, residual, residual);
		auto m2 = operations::overlap_diagonal(st, basis, phi, residual);
		auto m3 = operations::overlap_diagonal(st, basis, residual, hresidual);
		auto m4 = operations::overlap_diagonal(st, basis, phi, hresidual);


		//DATAOPERATIONS
		for(int ist = 0; ist < st.num_states(); ist++){
			double ca = real(m1[ist])*real(m4[ist]) - real(m3[ist])*real(m2[ist]);
			double cb = real(norm[ist])*real(m3[ist]) - real(eigenvalues[ist])*real(m1[ist]);
			double cc = real(eigenvalues[ist])*real(m2[ist]) - real(m4[ist])*real(norm[ist]);

			lambda[ist] = 2.0*cc/(cb + sqrt(cb*cb - 4.0*ca*cc));
		}

		operations::shift(st, basis, lambda, residual, phi);		
		
	}

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/plane_wave.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("solvers::steepest_descent", "[steepest_descent]") {


}


#endif


#endif
