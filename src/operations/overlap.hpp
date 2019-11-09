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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math/array.hpp>
#include <cassert>
#include <multi/adaptors/blas.hpp>
#include <operations/integral.hpp>

#ifdef HAVE_CUDA
#include <multi/memory/adaptors/cuda/allocator.hpp>
#endif

namespace operations {

	template <class field_set_type>
  auto overlap(const field_set_type & phi1, const field_set_type & phi2){

		math::array<typename field_set_type::value_type, 2>  overlap_matrix({phi1.set_size(), phi1.set_size()});

		boost::multi::blas::gemm('N', 'C', phi1.basis().volume_element(), phi1, phi2, 0.0, overlap_matrix);

		return overlap_matrix;		
  }

	template <class field_set_type>
	auto overlap(const field_set_type & phi){
		return overlap(phi, phi);
	}

	template <class field_set_type>
	math::array<typename field_set_type::element_type, 1> overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

		using type = typename field_set_type::element_type;
		
		math::array<type, 1> overlap_vector(phi1.set_size());

		assert(size(overlap_vector) == phi1.set_size());

		//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifndef HAVE_CUDA

		//OPTIMIZATION: this can be done more efficiently
    for(int ii = 0; ii < phi1.set_size(); ii++){
			type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().size(); ip++) aa += conj(phi1[ip][ii])*phi2[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
    }

#else

		{
			auto npoints = phi1.basis().size();
			auto nst = phi1.set_size();
			auto vol_element = phi1.basis().volume_element();
			auto phi1p = begin(phi1);
			auto phi2p = begin(phi2);
			auto overlap = begin(overlap_vector);
			
			//OPTIMIZATION: here we should parallelize over points as well 
			gpu::run(phi1.set_size(),
							 [=] __device__ (auto ist){
								 type aa = 0.0;
								 for(int ip = 0; ip < npoints; ip++){
									 auto p1 = phi1p[ip][ist];
									 auto p2 = phi2p[ip][ist];
									 aa += conj(p1)*p2;
									 
								 }
								 
								 overlap[ist] = vol_element*aa;
							 });
		}
		
#endif
		
		return overlap_vector;		
  }
	
	template <class field_set_type>
	auto overlap_diagonal(const field_set_type & phi){
		
		return overlap_diagonal(phi, phi);
	}
	
	template <class field_type>
	auto overlap_single(const field_type & phi1, const field_type & phi2){
		return integral(phi1, phi2, [](auto t1, auto t2){ return conj(t1)*t2; });
	}
	
	template <class field_type>
	auto overlap_single(field_type & phi){
		return overlap_single(phi, phi);
	}
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

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

			REQUIRE(real(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
			REQUIRE(imag(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*N*(N + 1.0)*bas.volume_element()));
		
		}
		
}


#endif
#endif
