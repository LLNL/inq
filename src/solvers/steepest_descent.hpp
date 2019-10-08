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
	void steepest_descent(const operator_type & ham, const preconditioner_type & prec, basis::field_set<basis::fourier_space, field_set_type> & phi){

		const int num_steps = 5;

		for(int istep = 0; istep < num_steps; istep++){
			
			//calculate the residual
			
			auto residual = ham(phi);
			
			auto eigenvalues = operations::overlap_diagonal(residual, phi);
			auto norm =	operations::overlap_diagonal(phi);
			
			auto lambda = eigenvalues;
			
			//DATAOPERATIONS STL
			std::transform(lambda.begin(), lambda.end(), norm.begin(), lambda.begin(), [](auto lam, auto nor){ return lam /= -nor; });

			operations::shift(lambda, phi, residual);
			
			prec(residual);
			
			//now calculate the step size
			auto hresidual = ham(residual);

			auto mm = boost::multi::array<field_set_type, 2>({6, phi.set_size()});

			mm[0] = operations::overlap_diagonal(residual, residual);
			mm[1] = operations::overlap_diagonal(phi, residual);
			mm[2] = operations::overlap_diagonal(residual, hresidual);
			mm[3] = operations::overlap_diagonal(phi, hresidual);
			mm[4] = eigenvalues;
			mm[5] = norm;
			
			//DATAOPERATIONS STL
			std::transform(mm.rotated().begin(), mm.rotated().end(), lambda.begin(),
										 [](auto mm){
											 auto ca = real(mm[0]*mm[3] - mm[2]*mm[1]);
											 auto cb = real(mm[5]*mm[2] - mm[4]*mm[0]);
											 auto cc = real(mm[4]*mm[1] - mm[3]*mm[5]);
											 
											 return 2.0*cc/(cb + sqrt(cb*cb - 4.0*ca*cc));
										 });
			
			operations::shift(lambda, residual, phi);		

		}

		operations::orthogonalization(phi);
		
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
