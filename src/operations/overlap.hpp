/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP
#define OPERATIONS__OVERLAP

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <multi/array.hpp>
#include <basis/field_set.hpp>
#include <cassert>
#include <multi/adaptors/blas.hpp>

#ifdef HAVE_CUDA
#include <multi/memory/adaptors/cuda/allocator.hpp>
#endif

namespace operations {

	template <class field_set_type>
  auto overlap(const field_set_type & phi1, const field_set_type & phi2){

		boost::multi::array<typename field_set_type::value_type, 2>  overlap_matrix({phi1.set_size(), phi1.set_size()});

		boost::multi::blas::gemm('N', 'C', phi1.basis().volume_element(), phi1, phi2, 0.0, overlap_matrix);

		return overlap_matrix;		
  }

	template <class field_set_type>
	auto overlap(const field_set_type & phi){
		return overlap(phi, phi);
	}

#ifdef HAVE_CUDA
	template <class type>
	__global__ void overlap_diagonal_kernel(const long npoints, const int nst, const double vol_element,
																					const type * phi1, const type * phi2, type * overlap){
		
		int ist = blockIdx.x*blockDim.x + threadIdx.x;

		type aa = 0.0;
		for(int ip = 0; ip < npoints; ip++){
			auto p1 = phi1[ip*nst + ist];
			auto p2 = phi2[ip*nst + ist];
			aa += conj(p1)*p2;
		}
		
		overlap[ist] = vol_element*aa;

	}
#endif
	
	template <class field_set_type>
  auto overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

		namespace multi = boost::multi;

		using value_type = typename field_set_type::value_type;
		
		multi::array<value_type, 1>  overlap_vector(phi1.set_size());

		assert(size(overlap_vector) == phi1.set_size());

#ifndef HAVE_CUDA

		//DATAOPERATIONS
		
		//OPTIMIZATION: this can be done more efficiently
    for(int ii = 0; ii < phi1.set_size(); ii++){
			typename field_set_type::value_type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().size(); ip++) aa += conj(phi1[ip][ii])*phi2[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
    }

#else

		namespace cuda = multi::memory::cuda;
		
		multi::array<value_type, 1, cuda::allocator<complex>> overlap_cuda(phi1.set_size());
		multi::array<value_type, 2, cuda::allocator<complex>> phi1_cuda = phi1;
		multi::array<value_type, 2, cuda::allocator<complex>> phi2_cuda = phi2;

		//OPTIMIZATION: here we should parallelize over points as well 
		overlap_diagonal_kernel<value_type><<<1, phi1.set_size()>>>(phi1.basis().size(), phi1.set_size(), phi1.basis().volume_element(),
																																static_cast<value_type const *>(phi1_cuda.data()),
																																static_cast<value_type const *>(phi2_cuda.data()),
																																static_cast<value_type *>(overlap_cuda.data()));
		
		overlap_vector = overlap_cuda;
		
#endif
		
		return overlap_vector;		
  }
	
	template <class field_set_type>
	auto overlap_diagonal(const field_set_type & phi){
		
		return overlap_diagonal(phi, phi);
	}
	
	template <class field_type>
	auto overlap_single(const field_type & phi1, const field_type & phi2){
		
		//DATAOPERATIONS
		//OPTIMIZATION: this can be done with BLAS
		typename field_type::value_type overlap = 0.0;
		for(int ip = 0; ip < phi1.basis().size(); ip++) overlap += conj(phi1[ip])*phi2[ip];
		return overlap*phi1.basis().volume_element();
	}
	
	template <class field_type>
	auto overlap_single(field_type & phi){
		return overlap_single(phi, phi);
	}
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("function operations::overlap", "[operations::overlap]") {

	using namespace Catch::literals;

		const int N = 100;
		const int M = 12;
			
		basis::trivial bas(N);

		SECTION("double"){
		
			basis::field_set<basis::trivial, double> aa(bas, M);
			basis::field_set<basis::trivial, double> bb(bas, M);

			for(int ii = 0; ii < N; ii++){
				for(int jj = 0; jj < M; jj++){
					aa[ii][jj] = 20.0*(ii + 1)*sqrt(jj);
					bb[ii][jj] = -0.05/(ii + 1)*sqrt(jj);
				}
			}

			{
				auto cc = operations::overlap(aa, bb);
				
				for(int ii = 0; ii < M; ii++){
					for(int jj = 0; jj < M; jj++) REQUIRE(cc[ii][jj] == Approx(-sqrt(jj)*sqrt(ii)));
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa, bb);
				
				for(int jj = 0; jj < M; jj++) REQUIRE(dd[jj] == Approx(-jj));
			}
			
			for(int ii = 0; ii < N; ii++){
				for(int jj = 0; jj < M; jj++){
					aa[ii][jj] = sqrt(ii)*sqrt(jj);
				}
			}

			{
				auto cc = operations::overlap(aa);
								
				for(int ii = 0; ii < M; ii++){
					for(int jj = 0; jj < M; jj++) REQUIRE(cc[ii][jj] == Approx(0.5*N*(N - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa);
								
				for(int jj = 0; jj < M; jj++) REQUIRE(dd[jj] == Approx(0.5*N*(N - 1.0)*bas.volume_element()*jj));
			}
					
			
		}

		SECTION("complex"){
		
			basis::field_set<basis::trivial, complex> aa(bas, M);
			basis::field_set<basis::trivial, complex> bb(bas, M);

			for(int ii = 0; ii < N; ii++){
				for(int jj = 0; jj < M; jj++){
					aa[ii][jj] = 20.0*(ii + 1)*sqrt(jj)*exp(complex(0.0, -M_PI/4 + M_PI/7*ii));
					bb[ii][jj] = -0.05/(ii + 1)*sqrt(jj)*exp(complex(0.0, M_PI/4 + M_PI/7*ii));
				}
			}

			{
				auto cc = operations::overlap(aa, bb);

				for(int ii = 0; ii < M; ii++){
					for(int jj = 0; jj < M; jj++) {
						REQUIRE(fabs(real(cc[ii][jj])) < 1.0e-14);
						REQUIRE(imag(cc[ii][jj]) == Approx(sqrt(jj)*sqrt(ii)));
					}
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa, bb);
				
				for(int jj = 0; jj < M; jj++){
					REQUIRE(fabs(real(dd[jj])) < 1.0e-14);
					REQUIRE(imag(dd[jj]) == Approx(-jj));
				}
			}
			
			for(int ii = 0; ii < N; ii++){
				for(int jj = 0; jj < M; jj++){
					aa[ii][jj] = sqrt(ii)*sqrt(jj)*exp(complex(0.0, M_PI/65.0*ii));
				}
			}

			{
				auto cc = operations::overlap(aa);

				std::cout << cc[0][0] << '\t' << cc[0][1] << '\t' << cc[0][2] << '\t' << cc[0][3]<< std::endl;
				std::cout << cc[1][0] << '\t' << cc[1][1] << '\t' << cc[1][2] << '\t' << cc[1][3]<< std::endl;
				std::cout << cc[2][0] << '\t' << cc[2][1] << '\t' << cc[2][2] << '\t' << cc[2][3]<< std::endl;
				std::cout << cc[3][0] << '\t' << cc[3][1] << '\t' << cc[3][2] << '\t' << cc[3][3]<< std::endl;
				
				for(int ii = 0; ii < M; ii++){
					for(int jj = 0; jj < M; jj++){
						REQUIRE(real(cc[ii][jj]) == Approx(0.5*N*(N - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
						REQUIRE(fabs(imag(cc[ii][jj])) < 1e-13);
					}
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa);
								
				for(int jj = 0; jj < M; jj++) REQUIRE(real(dd[jj]) == Approx(0.5*N*(N - 1.0)*bas.volume_element()*jj));
			}
					
			
		}

		
		SECTION("Overlap single double"){
			
			basis::field<basis::trivial, double> aa(bas);
			basis::field<basis::trivial, double> bb(bas);
			
			aa = 2.0;
			bb = 0.8;
		
			REQUIRE(operations::overlap_single(aa, bb) == 1.6_a);
			
			for(int ii = 0; ii < N; ii++)	{
				aa[ii] = pow(ii + 1, 2);
				bb[ii] = 1.0/(ii + 1);
			}
			
			REQUIRE(operations::overlap_single(aa, bb) == Approx(0.5*N*(N + 1.0)*bas.volume_element()));
			
		}
		
		SECTION("Integral product complex"){
			
			basis::field<basis::trivial, complex> aa(bas);
			basis::field<basis::trivial, complex> bb(bas);
			
			aa = complex(2.0, -0.3);
			bb = complex(0.8, 0.01);
		
			REQUIRE(real(operations::overlap_single(aa, bb)) == 1.597_a);
			REQUIRE(imag(operations::overlap_single(aa, bb)) == 0.26_a);
		
			for(int ii = 0; ii < N; ii++)	{
				aa[ii] = pow(ii + 1, 2)*exp(complex(0.0, -M_PI/8 + 2.0*M_PI/(ii + 1)));
				bb[ii] = 1.0/(ii + 1)*exp(complex(0.0, M_PI/8 + 2.0*M_PI/(ii + 1)));
			}

			std::cout << operations::overlap_single(aa, bb) << std::endl;
			
			REQUIRE(real(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
			REQUIRE(imag(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
		
	}

}


#endif
#endif
