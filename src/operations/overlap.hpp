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
#ifdef HAVE_CUDA
#include "multi/adaptors/blas/cuda.hpp" // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>
#include <operations/integral.hpp>

namespace operations {

	template <class field_set_type>
  auto overlap(const field_set_type & phi1, const field_set_type & phi2){

		using boost::multi::blas::gemm;
		using boost::multi::blas::hermitized;

		return gemm(phi1.basis().volume_element(), hermitized(phi2.matrix()), phi1.matrix());

  }

	template <class field_set_type>
	auto overlap(const field_set_type & phi){
		using boost::multi::blas::herk;
		using boost::multi::blas::hermitized;
		
		return herk(phi.basis().volume_element(), hermitized(phi.matrix()));

	}

	template <class field_set_type>
	math::array<typename field_set_type::element_type, 1> overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

		using type = typename field_set_type::element_type;
		
		math::array<type, 1> overlap_vector(phi1.set_dist().local_size());

		assert(size(overlap_vector) == phi1.set_dist().local_size());

		//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifndef HAVE_CUDA

		//OPTIMIZATION: this can be done more efficiently
    for(int ii = 0; ii < phi1.set_dist().local_size(); ii++){
			type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().dist().local_size(); ip++) aa += conj(phi1.matrix()[ip][ii])*phi2.matrix()[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
    }

#else

		{
			auto npoints = phi1.basis().dist().local_size();
			auto vol_element = phi1.basis().volume_element();
			auto phi1p = begin(phi1.matrix());
			auto phi2p = begin(phi2.matrix());
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

		if(phi1.basis().dist().parallel()){
			phi1.basis_comm().all_reduce_in_place_n(static_cast<type *>(overlap_vector.data()), overlap_vector.size(), std::plus<>{});
		}
		
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

		const int npoint = 100;
		const int nvec = 12;
			
		basis::trivial bas(npoint);

		SECTION("double"){
		
			basis::field_set<basis::trivial, double> aa(bas, nvec);
			basis::field_set<basis::trivial, double> bb(bas, nvec);

			for(int ii = 0; ii < npoint; ii++){
				for(int jj = 0; jj < nvec; jj++){
					aa.matrix()[ii][jj] = 20.0*(ii + 1)*sqrt(jj);
					bb.matrix()[ii][jj] = -0.05/(ii + 1)*sqrt(jj);
				}
			}

			{
				auto cc = operations::overlap(aa, bb);

				for(int ii = 0; ii < nvec; ii++){
					for(int jj = 0; jj < nvec; jj++) REQUIRE(cc[ii][jj] == Approx(-sqrt(jj)*sqrt(ii)));
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa, bb);
				
				for(int jj = 0; jj < nvec; jj++) REQUIRE(dd[jj] == Approx(-jj));
			}
			
			for(int ii = 0; ii < npoint; ii++){
				for(int jj = 0; jj < nvec; jj++){
					aa.matrix()[ii][jj] = sqrt(ii)*sqrt(jj);
				}
			}

			{
				auto cc = operations::overlap(aa);

				REQUIRE(std::get<0>(sizes(cc)) == nvec);
				REQUIRE(std::get<1>(sizes(cc)) == nvec);
				
				for(int ii = 0; ii < nvec; ii++){
					for(int jj = 0; jj < nvec; jj++) REQUIRE(cc[ii][jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa);
								
				for(int jj = 0; jj < nvec; jj++) REQUIRE(dd[jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
			}
					
			
		}

		SECTION("complex"){
		
			basis::field_set<basis::trivial, complex> aa(bas, nvec);
			basis::field_set<basis::trivial, complex> bb(bas, nvec);

			for(int ii = 0; ii < npoint; ii++){
				for(int jj = 0; jj < nvec; jj++){
					aa.matrix()[ii][jj] = 20.0*(ii + 1)*sqrt(jj)*exp(complex(0.0, -M_PI/4 + M_PI/7*ii));
					bb.matrix()[ii][jj] = -0.05/(ii + 1)*sqrt(jj)*exp(complex(0.0, M_PI/4 + M_PI/7*ii));
				}
			}

			{
				auto cc = operations::overlap(aa, bb);

				REQUIRE(std::get<0>(sizes(cc)) == nvec);
				REQUIRE(std::get<1>(sizes(cc)) == nvec);
				
				for(int ii = 0; ii < nvec; ii++){
					for(int jj = 0; jj < nvec; jj++) {
						REQUIRE(fabs(real(cc[ii][jj])) < 1.0e-14);
						REQUIRE(imag(cc[ii][jj]) == Approx(sqrt(jj)*sqrt(ii)));
					}
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa, bb);

				REQUIRE(std::get<0>(sizes(dd)) == nvec);
				
				for(int jj = 0; jj < nvec; jj++){
					REQUIRE(fabs(real(dd[jj])) < 1.0e-14);
					REQUIRE(imag(dd[jj]) == Approx(-jj));
				}
			}
			
			for(int ii = 0; ii < npoint; ii++){
				for(int jj = 0; jj < nvec; jj++){
					aa.matrix()[ii][jj] = sqrt(ii)*sqrt(jj)*exp(complex(0.0, M_PI/65.0*ii));
				}
			}

			{
				auto cc = operations::overlap(aa);

				REQUIRE(std::get<0>(sizes(cc)) == nvec);
				REQUIRE(std::get<1>(sizes(cc)) == nvec);
				
				for(int ii = 0; ii < nvec; ii++){
					for(int jj = 0; jj < nvec; jj++){
						REQUIRE(real(cc[ii][jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
						REQUIRE(fabs(imag(cc[ii][jj])) < 1e-13);
					}
				}
			}

			{
				auto dd = operations::overlap_diagonal(aa);

				REQUIRE(std::get<0>(sizes(dd)) == nvec);
				
				for(int jj = 0; jj < nvec; jj++) REQUIRE(real(dd[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
			}
					
			
		}

		
		SECTION("Overlap single double"){
			
			basis::field<basis::trivial, double> aa(bas);
			basis::field<basis::trivial, double> bb(bas);
			
			aa = 2.0;
			bb = 0.8;
		
			REQUIRE(operations::overlap_single(aa, bb) == 1.6_a);
			
			for(int ii = 0; ii < npoint; ii++)	{
				aa.linear()[ii] = pow(ii + 1, 2);
				bb.linear()[ii] = 1.0/(ii + 1);
			}
			
			REQUIRE(operations::overlap_single(aa, bb) == Approx(0.5*npoint*(npoint + 1.0)*bas.volume_element()));
			
		}
		
		SECTION("Integral product complex"){
			
			basis::field<basis::trivial, complex> aa(bas);
			basis::field<basis::trivial, complex> bb(bas);
			
			aa = complex(2.0, -0.3);
			bb = complex(0.8, 0.01);
		
			REQUIRE(real(operations::overlap_single(aa, bb)) == 1.597_a);
			REQUIRE(imag(operations::overlap_single(aa, bb)) == 0.26_a);
		
			for(int ii = 0; ii < npoint; ii++)	{
				aa.linear()[ii] = pow(ii + 1, 2)*exp(complex(0.0, -M_PI/8 + 2.0*M_PI/(ii + 1)));
				bb.linear()[ii] = 1.0/(ii + 1)*exp(complex(0.0, M_PI/8 + 2.0*M_PI/(ii + 1)));
			}

			REQUIRE(real(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*npoint*(npoint + 1.0)*bas.volume_element()));
			REQUIRE(imag(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*npoint*(npoint + 1.0)*bas.volume_element()));
		
		}

		SECTION("complex 1x1"){

			const int nvec = 1;
			const int npoint = 60000;
			
			basis::field_set<basis::trivial, complex> aa(bas, nvec);

			{
				auto cc = operations::overlap(aa);

				REQUIRE(std::get<0>(sizes(cc)) == nvec);
				REQUIRE(std::get<1>(sizes(cc)) == nvec);
				
				for(int ii = 0; ii < nvec; ii++){
					for(int jj = 0; jj < nvec; jj++){
						REQUIRE(real(cc[ii][jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
						REQUIRE(fabs(imag(cc[ii][jj])) < 1e-13);
					}
				}
			}

		}

}


#endif
#endif
