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
#include <math/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/shift.hpp>

namespace solvers {

	template <class operator_type, class preconditioner_type, class field_set_type>
	void steepest_descent(const operator_type & ham, const preconditioner_type & prec, field_set_type & phi){

		const int num_steps = 5;

		for(int istep = 0; istep < num_steps; istep++){
			
			//calculate the residual
			
			auto residual = ham(phi);
			
			auto eigenvalues = operations::overlap_diagonal(residual, phi);
			auto norm =	operations::overlap_diagonal(phi);
			
			auto lambda = eigenvalues;
			
			//DATAOPERATIONS STL + GPU::RUN TRANSFORM
#ifdef HAVE_CUDA
			gpu::run(phi.set_size(),
							 [lam = begin(lambda), nor = begin(norm)] __device__
							 (auto ist){
								 lam[ist] = lam[ist]/(-1.0*nor[ist]);
							 });
#else
			std::transform(lambda.begin(), lambda.end(), norm.begin(), lambda.begin(), [](auto lam, auto nor){ return lam /= -nor; });
#endif
			
			operations::shift(lambda, phi, residual);

			prec(residual);
			
			//now calculate the step size
			auto hresidual = ham(residual);

			auto mm = math::array<typename field_set_type::element_type, 2>({6, phi.set_size()});

			mm[0] = operations::overlap_diagonal(residual, residual);
			mm[1] = operations::overlap_diagonal(phi, residual);
			mm[2] = operations::overlap_diagonal(residual, hresidual);
			mm[3] = operations::overlap_diagonal(phi, hresidual);
			mm[4] = eigenvalues;
			mm[5] = norm;
			
			//DATAOPERATIONS STL + GPU::RUN TRANSFORM
#ifdef HAVE_CUDA
			gpu::run(phi.set_size(),
							 [m = begin(mm), lam = begin(lambda)]
							 __device__ (auto ist){
								 auto ca = real(m[0][ist]*m[3][ist] - m[2][ist]*m[1][ist]);
								 auto cb = real(m[5][ist]*m[2][ist] - m[4][ist]*m[0][ist]);
								 auto cc = real(m[4][ist]*m[1][ist] - m[3][ist]*m[5][ist]);

								 auto den = cb + sqrt(cb*cb - 4.0*ca*cc);

								 if(fabs(den) < 1e-15) lam[ist] = complex(0.0, 0.0); //this happens if we are perfectly converged
								 lam[ist] = complex(2.0*cc/den, 0.0);
							 });
#else
			std::transform(mm.rotated().begin(), mm.rotated().end(), lambda.begin(),
										 [](auto mm){
											 auto ca = real(mm[0]*mm[3] - mm[2]*mm[1]);
											 auto cb = real(mm[5]*mm[2] - mm[4]*mm[0]);
											 auto cc = real(mm[4]*mm[1] - mm[3]*mm[5]);

											 auto den = cb + sqrt(cb*cb - 4.0*ca*cc); 

											 if(fabs(den) < 1e-15) return 0.0;  //this happens if we are perfectly converged
											 return 2.0*cc/den;
										 });
#endif
			
			operations::shift(lambda, residual, phi);

		}

		operations::orthogonalization(phi);
		
	}

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("solvers::steepest_descent", "[solvers::steepest_descent]") {

  const int npoint = 100;
  const int nvec = 12;
  
  basis::trivial bas(npoint);
	
	math::array<complex, 2> identity_matrix({npoint, npoint});
  
	for(int ip = 0; ip < npoint; ip++){
		for(int jp = 0; jp < npoint; jp++){
			identity_matrix[ip][jp] = 0.0;
			if(ip == jp) identity_matrix[ip][jp] = 1.0;
		}
	}
  
	operations::matrix_operator<math::array<complex, 2>> identity(std::move(identity_matrix));

	SECTION("Diagonal matrix complex"){
  
    math::array<complex, 2> diagonal_matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        diagonal_matrix[ip][jp] = 0.0;
        if(ip == jp) diagonal_matrix[ip][jp] = ip;
      }
    }
    
    operations::matrix_operator<math::array<complex, 2>> diagonal_op(std::move(diagonal_matrix));
    
    basis::field_set<basis::trivial, complex> phi(bas, nvec);

		for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        phi.matrix()[ip][ivec] = exp(complex(0.0, (ip*ivec)*0.1));
      }
    }

		operations::orthogonalization(phi);
		
		for(int iter = 0; iter < 100; iter++){

			tfm::format(std::cout, "  Iteration %4d:\n", iter);
			
			solvers::steepest_descent(diagonal_op, identity, phi);
			
			auto residual = diagonal_op(phi);
			auto eigenvalues = operations::overlap_diagonal(phi, residual);
			operations::shift(eigenvalues, phi, residual, -1.0);
			auto normres = operations::overlap_diagonal(residual);
			
			for(int ivec = 0; ivec < phi.set_size(); ivec++){
				tfm::format(std::cout, "    state %4d  evalue = %18.12f  res = %15.10e\n", ivec + 1, real(eigenvalues[ivec]), real(normres[ivec]));
			}
		}
 	
 }

#if 0
	SECTION("Periodic Laplacian matrix complex"){
  
    math::array<complex, 2> laplacian_matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        laplacian_matrix[ip][jp] = 0.0;
				identity_matrix[ip][jp] = 0.0;
        if(ip == jp) laplacian_matrix[ip][jp] = -1.0;
        if(ip == jp + 1 or ip == jp - 1) laplacian_matrix[ip][jp] = 2.0;
      }
    }
    //the periodic part
    laplacian_matrix[0][npoint - 1] = 2.0/bas.volume_element();
    laplacian_matrix[npoint - 1][0] = 2.0/bas.volume_element();
    
    operations::matrix_operator<math::array<complex, 2>> laplacian(std::move(laplacian_matrix));
    
    basis::field_set<basis::trivial, complex> phi(bas, nvec);

		phi = 0.0;
		for(int ivec = 0; ivec < nvec; ivec++) phi.matrix()[ivec][ivec] = 1.0;

		for(int iter = 0; iter < 100; iter++){
			
			solvers::steepest_descent(laplacian, identity, phi);
			
			auto residual = laplacian(phi);
			auto eigenvalues = operations::overlap_diagonal(phi, residual);
			operations::shift(eigenvalues, phi, residual, -1.0);
			auto normres = operations::overlap_diagonal(residual);
			
			for(int ivec = 0; ivec < phi.set_size(); ivec++){
				tfm::format(std::cout, " state %4d  evalue = %18.12f  res = %5.0e\n", ivec + 1, real(eigenvalues[ivec]), real(normres[ivec]));
			}
		}
 	
 }
#endif
	

}


#endif


#endif
