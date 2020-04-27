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
#include <math/vec3d.hpp>
#include <math/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/shift.hpp>
#include <operations/orthogonalize.hpp>

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
			
			//DATAOPERATIONS GPU::RUN TRANSFORM
			gpu::run(phi.set_size(),
							 [lam = begin(lambda), nor = begin(norm)]
							 GPU_LAMBDA (auto ist){
								 lam[ist] = lam[ist]/(-real(nor[ist]));
							 });
			
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
			
			//DATAOPERATIONS GPU::RUN TRANSFORM
			gpu::run(phi.set_size(),
							 [m = begin(mm), lam = begin(lambda)]
							 GPU_LAMBDA (auto ist){
								 auto ca = real(m[0][ist]*m[3][ist] - m[2][ist]*m[1][ist]);
								 auto cb = real(m[5][ist]*m[2][ist] - m[4][ist]*m[0][ist]);
								 auto cc = real(m[4][ist]*m[1][ist] - m[3][ist]*m[5][ist]);

								 auto den = cb + sqrt(cb*cb - 4.0*ca*cc);

								 if(fabs(den) < 1e-15) lam[ist] = complex(0.0, 0.0); //this happens if we are perfectly converged
								 lam[ist] = complex(2.0*cc/den, 0.0);
							 });

			operations::shift(lambda, residual, phi);

		}

		operations::orthogonalize(phi);
		
	}

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>
#include <operations/matrix_operator.hpp>

TEST_CASE("solvers::steepest_descent", "[solvers::steepest_descent]") {

	using namespace Catch::literals;
	
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
  
	operations::matrix_operator<complex> identity(std::move(identity_matrix));

	SECTION("Diagonal matrix complex"){
  
    math::array<complex, 2> diagonal_matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        diagonal_matrix[ip][jp] = 0.0;
        if(ip == jp) diagonal_matrix[ip][jp] = ip + 1.0;
      }
    }
    
    operations::matrix_operator<complex> diagonal_op(std::move(diagonal_matrix));
    
    basis::field_set<basis::trivial, complex> phi(bas, nvec);

		for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        phi.matrix()[ip][ivec] = exp(complex(0.0, (ip*ivec)*0.1));
      }
    }

		operations::orthogonalize(phi);

		const int num_iter = 100;
		
		for(int iter = 0; iter < num_iter; iter++){

			solvers::steepest_descent(diagonal_op, identity, phi);
			
			auto residual = diagonal_op(phi);
			auto eigenvalues = operations::overlap_diagonal(phi, residual);
			operations::shift(eigenvalues, phi, residual, -1.0);
			auto normres = operations::overlap_diagonal(residual);
			
			/*
				tfm::format(std::cout, "  Iteration %4d:\n", iter);
				
				for(int ivec = 0; ivec < phi.set_size(); ivec++){
				tfm::format(std::cout, "    state %4d  evalue = %18.12f  res = %15.10e\n", ivec + 1, real(eigenvalues[ivec]), real(normres[ivec]));
				}
			*/

			if(num_iter - 1 == iter){

				CHECK(fabs(eigenvalues[0]) == 1.000000001634_a);
				CHECK(fabs(eigenvalues[1]) == 2.000000003689_a);
				CHECK(fabs(eigenvalues[2]) == 3.000000001955_a);
				CHECK(fabs(eigenvalues[3]) == 4.000000001760_a);
				CHECK(fabs(eigenvalues[4]) == 5.000000002561_a);
				CHECK(fabs(eigenvalues[5]) == 6.000000003127_a);
				CHECK(fabs(eigenvalues[6]) == 7.000000002312_a);
				CHECK(fabs(eigenvalues[7]) == 8.000000000292_a);
				CHECK(fabs(eigenvalues[8]) == 8.999999999033_a);
				CHECK(fabs(eigenvalues[9]) == 9.999999998497_a);
				CHECK(fabs(eigenvalues[10]) == 10.999999998768_a);
				CHECK(fabs(eigenvalues[11]) == 11.999999998422_a);

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
		
    operations::matrix_operator<complex> laplacian(std::move(laplacian_matrix));
    
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
