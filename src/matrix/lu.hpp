/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATRIX__LU
#define INQ__MATRIX__LU

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <matrix/gather_scatter.hpp>


#ifdef ENABLE_CUDA
#include <cusolverDn.h>
#endif

#include "FC.h"

#define dgetrf FC_GLOBAL(dgetrf, DGETRF) 
extern "C" void dgetrf(const int * m, const int * n, double * a, const int * lda, int * ipiv, int * info);

namespace inq {
namespace matrix {

template <class matrix_type>
void lu_raw(matrix_type && matrix){
  
	const int nst = matrix.size();
	int info;
  
#ifdef NO_ENABLE_CUDA
	{

		CALI_CXX_MARK_SCOPE("cuda_dgetrf");
		
		cusolverDnHandle_t cusolver_handle;
			
		[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
		assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
		//query the work size
		int lwork;
		cusolver_status = cusolverDnDgetrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(matrix.data_elements()), nst, &lwork);
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

		cusolver_status = cusolverDnDgetrf(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(matrix.data_elements()), nst, work, lwork, devInfo);
		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		cudaDeviceSynchronize();
		info = *devInfo ;
		
		cudaFree(work);
		cudaFree(devInfo);
		cusolverDnDestroy(cusolver_handle);
			
	}
#else
	{
		CALI_CXX_MARK_SCOPE("cuda_dgetrf");
		auto ipiv = (int *) malloc(nst*sizeof(int));
		dgetrf(&nst, &nst, raw_pointer_cast(matrix.data_elements()), &nst, ipiv, &info);
		free(ipiv);
															 
	}
#endif
	if(info != 0){
		std::printf("Error: Failed LU decomposition, info is %10i.\n", info);
		abort();
	}
  
}

template <typename DistributedMatrix>
void lu(DistributedMatrix & matrix) {

  assert(matrix.sizex() == matrix.sizey());
  
  auto full_matrix = matrix::gather(matrix, /* root = */ 0);
  if(matrix.comm().root()) lu_raw(full_matrix);
  matrix::scatter(full_matrix, matrix, /* root = */ 0);

}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_MATRIX_LU_UNIT_TEST
#undef INQ_MATRIX_LU_UNIT_TEST

#include <gpu/array.hpp>
#include <matrix/distributed.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
  using namespace Catch::literals;
  
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	parallel::cartesian_communicator<2> cart_comm(comm, {});


  SECTION("Double 2x2"){
    
		using namespace inq;
    using namespace Catch::literals;
		
		gpu::array<double, 2> array({2, 2});
		
		array[0][0] = 6432.12;
		array[0][1] = 4502.48;
		array[1][0] = 4502.48;
		array[1][1] = 3151.74;

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		matrix::lu(matrix);

    array = matrix::all_gather(matrix);

		//comparison with octave values, the results here are trasposed because lapack uses Fortran
		CHECK(real(array[0][0]) == 6.4321e+03_a);
		CHECK(real(array[1][0]) == 4.5025e+03_a);
		CHECK(real(array[0][1]) == 7.0000e-01_a);
		CHECK(real(array[1][1]) == 6.8000e-03_a);    
  }

}
#endif
