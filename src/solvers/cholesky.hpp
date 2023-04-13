/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__CHOLESKY
#define INQ__SOLVERS__CHOLESKY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/complex.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

#ifdef ENABLE_CUDA
#include <cusolverDn.h>
#endif

#include "FC.h"

#define zpotrf FC_GLOBAL(zpotrf, ZPOTRF) 
extern "C" void zpotrf(const char * uplo, const int * n, inq::complex * a, const int * lda, int * info);

namespace inq {
namespace solvers {

template <class matrix_type>
void cholesky(matrix_type && matrix, bool nocheck = false){
  
	const int nst = matrix.size();
	int info;
  
#ifdef ENABLE_CUDA
	{

		CALI_CXX_MARK_SCOPE("cuda_zpotrf");
		
		cusolverDnHandle_t cusolver_handle;
			
		[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
		assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
		//query the work size
		int lwork;
		cusolver_status = cusolverDnZpotrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(matrix.data_elements()), nst, &lwork);
		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		assert(lwork >= 0);
			
		//allocate the work array
		cuDoubleComplex * work;
		[[maybe_unused]] auto cuda_status = cudaMalloc((void**)&work, sizeof(cuDoubleComplex)*lwork);
		assert(cudaSuccess == cuda_status);

		//finaly do the decomposition
		int * devInfo;
		cuda_status = cudaMallocManaged((void**)&devInfo, sizeof(int));
		assert(cudaSuccess == cuda_status);

		cusolver_status = cusolverDnZpotrf(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(matrix.data_elements()), nst, work, lwork, devInfo);
		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		cudaDeviceSynchronize();
		info = *devInfo ;
		
		cudaFree(work);
		cudaFree(devInfo);
		cusolverDnDestroy(cusolver_handle);
			
	}
#else
	{
		CALI_CXX_MARK_SCOPE("cuda_zpotrf");
		zpotrf("U", &nst, raw_pointer_cast(matrix.data_elements()), &nst, &info);
	}
#endif

	if(not nocheck and info < 0){
		std::printf("Error: Failed orthogonalization in ZPOTRF! info is %10i.\n", info);
		abort();
	} else if(info > 0) {
		std::printf("Warning: Imperfect orthogonalization in ZPOTRF! info is %10i, subspace size is %10i\n", info, nst); 
	}

  gpu::run(nst,
           [mat = begin(matrix), nst] GPU_LAMBDA (auto ist){
             for(auto jst = ist + 1; jst < (decltype(jst)) nst; jst++) mat[ist][jst] = 0.0;
           });
  
}

}
}
#endif

#ifdef INQ_SOLVERS_CHOLESKY_UNIT_TEST
#undef INQ_SOLVERS_CHOLESKY_UNIT_TEST

#include <math/array.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	SECTION("Complex 2x2"){
    
		using namespace inq;
    using namespace Catch::literals;
		
		math::array<complex, 2> matrix({2, 2});
		
		matrix[0][0] = 6432.12;
		matrix[0][1] = 4502.48;
		matrix[1][0] = 4502.48;
		matrix[1][1] = 3151.74;
    
		solvers::cholesky(matrix);

		CHECK(real(matrix[0][0]) == 80.2005_a);
		CHECK(real(matrix[0][1]) == 0.0_a);
    CHECK(real(matrix[1][0]) == 56.1402992511_a);
    CHECK(real(matrix[1][1]) == 0.0824620974_a);    
		
  }
}
#endif
