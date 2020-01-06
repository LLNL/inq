/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__ORTHOGONALIZATION
#define OPERATIONS__ORTHOGONALIZATION

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa

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
#include <math/complex.hpp>
#include <basis/field_set.hpp>
#include <cstdlib>
#include <multi/adaptors/blas/trsm.hpp>

#ifdef HAVE_CUDA
#include <cusolverDn.h>
#endif

#define zpotrf FC_FUNC(zpotrf, ZPOTRF) 
extern "C" void zpotrf(const char * uplo, const int * n, complex * a, const int * lda, int * info);

//#define blas_ztrsm FC_FUNC(ztrsm, ZTRSM) 
//extern "C" void blas_ztrsm(const char& side, const char& uplo, const char& transa, const char& diag,
//											const long& m, const long& n, const complex& alpha, const complex * a, const long& lda, complex * B, const long& ldb);


namespace operations {

	template <class field_set_type>
  void orthogonalization(field_set_type & phi){

		auto olap = overlap(phi);

		const int nst = phi.set_size();
		
		//DATAOPERATIONS RAWLAPACK zpotrf
#ifdef HAVE_CUDA
		{
			cusolverDnHandle_t cusolver_handle;
			
			auto cusolver_status = cusolverDnCreate(&cusolver_handle);
			assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
			//query the work size
			int lwork;
			cusolver_status = cusolverDnZpotrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(olap.data()), nst, &lwork);
			assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
			assert(lwork >= 0);
			
			//allocate the work array
			cuDoubleComplex * work;
			auto cuda_status = cudaMalloc((void**)&work, sizeof(cuDoubleComplex)*lwork);
			assert(cudaSuccess == cuda_status);

			//finaly do the decomposition
			int * devInfo;
			cuda_status = cudaMallocManaged((void**)&devInfo, sizeof(int));
			assert(cudaSuccess == cuda_status);

			cusolver_status = cusolverDnZpotrf(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(olap.data()), nst, work, lwork, devInfo);
			assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
			cudaDeviceSynchronize();
			assert(*devInfo == 0);

			cudaFree(work);
			cudaFree(devInfo);
			cusolverDnDestroy(cusolver_handle);
			
		}
#else
		int info;
		zpotrf("U", &nst, olap.data(), &nst, &info);
		assert(info == 0);
#endif

		//DATAOPERATIONS trsm
		using boost::multi::blas::hermitized;
		using boost::multi::blas::filling;
		
		trsm(filling::lower, olap, hermitized(phi.matrix()));

  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::orthogonalization", "[operations::orthogonalization]") {

	using namespace Catch::literals;
	using math::d3vector;

	double ecut = 25.0;
	double ll = 6.3;

	ions::geometry geo;
	ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
	basis::real_space pw(cell, input::basis::cutoff_energy(ecut));

	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms());
	
	SECTION("Dimension 3"){
		basis::field_set<basis::real_space, complex> phi(pw, 3);
		
		operations::randomize(phi);
		
		operations::orthogonalization(phi);
		
		auto olap = operations::overlap(phi);
		
		std::cout << "------" << std::endl;
		
		std::cout << olap[0][0] << '\t' << olap[0][1] << '\t' << olap[0][2] << std::endl;
		std::cout << olap[1][0] << '\t' << olap[1][1] << '\t' << olap[1][2] << std::endl;
		std::cout << olap[2][0] << '\t' << olap[2][1] << '\t' << olap[2][2] << std::endl;
		
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					REQUIRE(real(olap[ii][ii]) == 1.0_a);
					REQUIRE(fabs(imag(olap[ii][ii])) < 1e-14);
			} else {
					REQUIRE(fabs(olap[ii][jj]) < 1e-14);
				}
			}
		}
	}

	SECTION("Dimension 100"){
		basis::field_set<basis::real_space, complex> phi(pw, 100);
		
		operations::randomize(phi);
		
		operations::orthogonalization(phi);
		
		auto olap = operations::overlap(phi);
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					REQUIRE(real(olap[ii][ii]) == 1.0_a);
					REQUIRE(fabs(imag(olap[ii][ii])) < 1e-14);
				} else {
					REQUIRE(fabs(olap[ii][jj]) < 1e-13);
				}
			}
		}
	}


	SECTION("Dimension 37 - double orthogonalization"){
		basis::field_set<basis::real_space, complex> phi(pw, 37);
		
		operations::randomize(phi);
		
		operations::orthogonalization(phi);
		operations::orthogonalization(phi);
		
		auto olap = operations::overlap(phi);
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					REQUIRE(real(olap[ii][ii]) == 1.0_a);
					REQUIRE(fabs(imag(olap[ii][ii])) < 1e-16);
				} else {
					REQUIRE(fabs(olap[ii][jj]) < 5e-16);
				}
			}
		}
	}

	
	
}


#endif

#endif
