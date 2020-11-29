/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__QRFACTORIZE
#define INQ__OPERATIONS__QRFACTORIZE

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

#include <inq_config.h>
#include <math/complex.hpp>
#include <basis/field_set.hpp>
#include <cstdlib>
#include <multi/adaptors/blas/trsm.hpp>

#ifdef ENABLE_CUDA
#include <cusolverDn.h>
#endif

#include "FC.h"

#define zgeqrf FC_GLOBAL(zgeqrf, ZGEQRF) 
extern "C" void zgeqrf(const int & m, const int & n,   inq::complex * a, const int & lda, inq::complex * tau, inq::complex * work, const int & lwork, int & info);
#define zungqr FC_GLOBAL(zungqr, ZUNGQR) 
extern "C" void zungqr(const int & m, const int & n, const int &k,  inq::complex * a, const int & lda, inq::complex * tau, inq::complex * work, const int & lwork, int & info);

namespace inq {
  namespace operations {

    template <class field_set_type>
    void qrfactorize(field_set_type & phi){

      assert(phi.set_comm().size() == 1);
      
      const int npw = std::get<0>(sizes(phi.matrix()));
      const int nst = phi.set_size();
      assert(npw >= nst);
      
      //DATAOPERATIONS RAWLAPACK zpotrf
#ifdef ENABLE_CUDA
      {
	std::printf("QR not implemented for CUDA!!!");
	exit(1);
	// 	cusolverDnHandle_t cusolver_handle;
			
	// 	[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
	// 	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
	// 	//query the work size
	// 	int lwork;
	// 	cusolver_status = cusolverDnZpotrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(olap.data()), nst, &lwork);
	// 	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	// 	assert(lwork >= 0);
			
	// 	//allocate the work array
	// 	cuDoubleComplex * work;
	// 	[[maybe_unused]] auto cuda_status = cudaMalloc((void**)&work, sizeof(cuDoubleComplex)*lwork);
	// 	assert(cudaSuccess == cuda_status);

	// 	//finaly do the decomposition
	// 	int * devInfo;
	// 	cuda_status = cudaMallocManaged((void**)&devInfo, sizeof(int));
	// 	assert(cudaSuccess == cuda_status);

	// 	cusolver_status = cusolverDnZpotrf(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(olap.data()), nst, work, lwork, devInfo);
	// 	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	// 	cudaDeviceSynchronize();
	// 	assert(*devInfo == 0);

	// 	cudaFree(work);
	// 	cudaFree(devInfo);
	// 	cusolverDnDestroy(cusolver_handle);		
      }
#else
      auto tphi = (complex *) malloc(npw*nst*sizeof(complex));
      
      //std::printf("before QR:\n");
      int count=0;
      for(int i=0; i<nst; i++){
	for(int j=0; j<npw; j++){
	  tphi[count]=phi.matrix()[j][i];
	  //std::printf("%16.8f ",tphi[count]);
	  count++ ;
	}
	//std::printf("\n");
      }
      int info;
      //math::array<complex, 1> tau(std::min(npw,nst));
      auto tau = (complex *) malloc(nst*sizeof(complex));
      inq::complex lwork_query;
      int lwork_q = -1;
      zgeqrf(npw, nst, tphi, npw, tau, &lwork_query, lwork_q, info);

      int lwork = int(real(lwork_query));
      //math::array<complex, 1> work(lwork);
      auto work = (complex *) malloc(lwork*sizeof(complex));
      for(int i=0;i<lwork;i++){
	work[i]=0.0;
      }
      zgeqrf(npw, nst, tphi, npw, tau, work, lwork, info);
      assert(info == 0);
      // for(int i=0;i<nst;i++){
      // 	for(int j=i;j<nst;j++){
      // 	  if(i==j) phi.matrix()[i][j]=1.0;
      // 	  else phi.matrix()[i][j]=0.0;
      // 	}
      // }
      lwork_q=-1;
      zungqr(npw, nst, nst, tphi, npw, tau, &lwork_query, lwork_q, info);
      lwork = int(real(lwork_query));
      //math::array<complex, 1> work2(lwork);
      auto work2 = (complex *) malloc(lwork*sizeof(complex));
      for(int i=0;i<lwork;i++){
	work2[i]=0.0;
      }
      zungqr(npw, nst, nst, tphi, npw, tau, work2, lwork, info);
      free(work2);
      free(work);
      free(tau);
      //std::printf("after zungqr:\n");
      count=0;
      for(int i=0; i<nst; i++){
	for(int j=0; j<npw; j++){
	  phi.matrix()[j][i]=tphi[count]/sqrt(phi.basis().volume_element());
	  //std::printf("%16.8f ",tphi[count]);
	  count++ ;
	}
	//std::printf("\n");
      }
      free(tphi);
      
      
#endif

    }
	

  }
}

#ifdef INQ_OPERATIONS_QRFACTORIZE_UNIT_TEST
#undef INQ_OPERATIONS_QRFACTORIZE_UNIT_TEST

#include <catch2/catch.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::qrfactorizeze", "[operations::qrfactorize]") {

	using namespace inq;
	using namespace Catch::literals;
	using math::vector3;

	double ecut = 25.0;
	double ll = 6.3;

	ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
	basis::real_space pw(cell, input::basis::cutoff_energy(ecut));

	SECTION("Dimension 3"){
		basis::field_set<basis::real_space, complex> phi(pw, 3);
		
		// operations::randomize(phi);
		
		// operations::orthogonalize(phi);
		
		// auto olap = operations::overlap(phi);
		
		// std::cout << "------" << std::endl;
		
		// std::cout << olap[0][0] << '\t' << olap[0][1] << '\t' << olap[0][2] << std::endl;
		// std::cout << olap[1][0] << '\t' << olap[1][1] << '\t' << olap[1][2] << std::endl;
		// std::cout << olap[2][0] << '\t' << olap[2][1] << '\t' << olap[2][2] << std::endl;
		
		
		// for(int ii = 0; ii < phi.set_size(); ii++){
		// 	for(int jj = 0; jj < phi.set_size(); jj++){
		// 		if(ii == jj) {
		// 			CHECK(real(olap[ii][ii]) == 1.0_a);
		// 			CHECK(fabs(imag(olap[ii][ii])) < 1e-14);
		// 	} else {
		// 			CHECK(fabs(olap[ii][jj]) < 1e-14);
		// 		}
		// 	}
		// }
	}
}


#endif

#endif
