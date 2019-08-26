/* -*- indent-tabs-mode: t -*- */

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

	template <class operator_type, class preconditioner_type, class field_set_type>
	void steepest_descent(const operator_type & ham, const preconditioner_type & prec, basis::field_set<basis::real_space, field_set_type> & phi){

		const int num_steps = 5;

		auto fphi = operations::space::to_fourier(phi);
		
		for(int istep = 0; istep < num_steps; istep++){
			
			//calculate the residual
			
			auto residual = ham(fphi);
			
			auto eigenvalues = operations::overlap_diagonal(residual, fphi);
			auto norm =	operations::overlap_diagonal(fphi);
			
			auto lambda(eigenvalues);
			
			//DATAOPERATIONS
			for(int ist = 0; ist < fphi.set_size(); ist++) lambda[ist] /= -norm[ist];
			
			operations::shift(lambda, fphi, residual);
			
			//OPTIMIZATIONS: precondition the residual here
			prec(residual);
			
			//now calculate the step size
			auto hresidual = ham(residual);
			
			auto m1 = operations::overlap_diagonal(residual, residual);
			auto m2 = operations::overlap_diagonal(fphi, residual);
			auto m3 = operations::overlap_diagonal(residual, hresidual);
			auto m4 = operations::overlap_diagonal(fphi, hresidual);
			
			
			//DATAOPERATIONS
			for(int ist = 0; ist < fphi.set_size(); ist++){
				double ca = real(m1[ist])*real(m4[ist]) - real(m3[ist])*real(m2[ist]);
				double cb = real(norm[ist])*real(m3[ist]) - real(eigenvalues[ist])*real(m1[ist]);
				double cc = real(eigenvalues[ist])*real(m2[ist]) - real(m4[ist])*real(norm[ist]);
				
				lambda[ist] = 2.0*cc/(cb + sqrt(cb*cb - 4.0*ca*cc));
			}
			
			operations::shift(lambda, residual, fphi);		

		}

		operations::orthogonalization(fphi);
		
		phi = operations::space::to_real(fphi);
		
	}

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("solvers::steepest_descent", "[steepest_descent]") {


}


#endif


#endif
