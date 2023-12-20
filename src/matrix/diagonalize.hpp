/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATRIX__DIAGONALIZE
#define INQ__MATRIX__DIAGONALIZE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <math/complex.hpp>
#include <matrix/gather_scatter.hpp>

#include <inq_config.h>

#ifdef ENABLE_CUDA
#include <cusolverDn.h>
#endif

#include "FC.h"

#include <utils/profiling.hpp>

#define dsyev FC_GLOBAL(dsyev, DZYEV)
extern "C" void dsyev(const char * jobz, const char * uplo, const int & n, double * a, const int & lda, double * w, double * work, const int & lwork, int & info);

#define zheev FC_GLOBAL(zheev, ZHEEV) 
extern "C" void zheev(const char * jobz, const char * uplo, const int & n, inq::complex * a, const int & lda, double * w, inq::complex * work, const int & lwork, double * rwork, int & info);

namespace inq {
namespace matrix {

template<class Alloc>
auto diagonalize_raw(gpu::array<double, 2, Alloc>& matrix){

	CALI_CXX_MARK_FUNCTION;

	// the matrix must be square
	assert(std::get<0>(sizes(matrix)) == std::get<1>(sizes(matrix)));

	int nn = std::get<0>(sizes(matrix));
    
	gpu::array<double, 1> eigenvalues(nn);

#ifdef ENABLE_CUDA
	{
		cusolverDnHandle_t cusolver_handle;
			
		[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
		assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
		//query the work size
		int lwork;
		cusolver_status = cusolverDnDsyevd_bufferSize(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, nn,
																									raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), &lwork);
		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		assert(lwork >= 0);

		//allocate the work array
		double * work;
		[[maybe_unused]] auto cuda_status = cudaMalloc((void**)&work, sizeof(double)*lwork);
		assert(cudaSuccess == cuda_status);

		//finally, diagonalize
		int * devInfo;
		cuda_status = cudaMallocManaged((void**)&devInfo, sizeof(int));
		assert(cudaSuccess == cuda_status);
			
		cusolver_status = cusolverDnDsyevd(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, nn,
																			 raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), work, lwork, devInfo);

		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		cudaDeviceSynchronize();
		assert(*devInfo == 0);

		cudaFree(work);
		cudaFree(devInfo);
		cusolverDnDestroy(cusolver_handle);
	}
#else
	double lwork_query;

	int info;
	dsyev("V", "U", nn, raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), &lwork_query, -1, info);

	int lwork = int(lwork_query);
	auto work = (double *) malloc(lwork*sizeof(complex));
    
	dsyev("V", "U", nn, raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), work, lwork, info);
	assert(info == 0);
		
	free(work);
#endif

	return eigenvalues;
}

template<class Alloc>
auto diagonalize_raw(gpu::array<complex, 2, Alloc>& matrix){

	CALI_CXX_MARK_FUNCTION;
	
	// the matrix must be square
	assert(std::get<0>(sizes(matrix)) == std::get<1>(sizes(matrix)));

	int nn = std::get<0>(sizes(matrix));
    
	gpu::array<double, 1> eigenvalues(nn);

#ifdef ENABLE_CUDA
	{
		cusolverDnHandle_t cusolver_handle;
			
		[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
		assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
		//query the work size
		int lwork;
		cusolver_status = cusolverDnZheevd_bufferSize(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, nn,
																									(cuDoubleComplex const *) raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), &lwork);
		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		assert(lwork >= 0);

		//allocate the work array
		cuDoubleComplex * work;
		[[maybe_unused]] auto cuda_status = cudaMalloc((void**)&work, sizeof(cuDoubleComplex)*lwork);
		assert(cudaSuccess == cuda_status);

		//finally, diagonalize
		int * devInfo;
		cuda_status = cudaMallocManaged((void**)&devInfo, sizeof(int));
		assert(cudaSuccess == cuda_status);
			
		cusolver_status = cusolverDnZheevd(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, nn,
																			 (cuDoubleComplex *) raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), work, lwork, devInfo);
		assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
		cudaDeviceSynchronize();
		assert(*devInfo == 0);
			
		cudaFree(devInfo);
		cusolverDnDestroy(cusolver_handle);
	}
#else
	complex lwork_query;

	auto rwork = (double *) malloc(std::max(1, 3*nn - 2)*sizeof(double));
    
	int info;
	zheev("V", "U", nn, raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), &lwork_query, -1, rwork, info);


	int lwork = int(real(lwork_query));
	auto work = (complex *) malloc(lwork*sizeof(complex));
    
	zheev("V", "U", nn, raw_pointer_cast(matrix.data_elements()), nn, raw_pointer_cast(eigenvalues.data_elements()), work, lwork, rwork, info);
    
	free(rwork);
	free(work);
#endif
		
	return eigenvalues;
}

template <typename DistributedMatrix>
auto diagonalize(DistributedMatrix & matrix) {

  assert(matrix.sizex() == matrix.sizey());
  
  gpu::array<double, 1> eigenvalues;

  auto full_matrix = matrix::gather(matrix, /* root = */ 0);
  
  if(matrix.comm().root()) {
    eigenvalues = diagonalize_raw(full_matrix);
  } else {
    eigenvalues = gpu::array<double, 1>(matrix.sizex());
  }
  
  assert(eigenvalues.size() == matrix.sizex());
  
	matrix::scatter(full_matrix, matrix, /* root = */ 0);
  
  matrix.comm().broadcast_n(raw_pointer_cast(eigenvalues.data_elements()), eigenvalues.num_elements(), 0);
  matrix.comm().barrier();
  
  return eigenvalues;
}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_MATRIX_DIAGONALIZE_UNIT_TEST
#undef INQ_MATRIX_DIAGONALIZE_UNIT_TEST

#include <gpu/array.hpp>
#include <matrix/distributed.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
  using namespace Catch::literals;
  
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	parallel::cartesian_communicator<2> cart_comm(comm, {});
	
	SECTION("Real diagonal 2x2"){
		
		gpu::array<double, 2> array({2, 2});
		
		array[0][0] = 4.0;
		array[0][1] = 0.0;
		array[1][0] = 0.0;
		array[1][1] = 2.0;

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		auto evalues = matrix::diagonalize(matrix);

    array = matrix::all_gather(matrix);

		CHECK(array[0][0] == 0.0_a);
		CHECK(array[0][1] == 1.0_a);
		CHECK(array[1][0] == 1.0_a);
		CHECK(array[0][0] == 0.0_a);

    CHECK(evalues[0] == 2.0_a);
    CHECK(evalues[1] == 4.0_a);

	}
	
	SECTION("Complex diagonal 2x2"){
		
		gpu::array<complex, 2> array({2, 2});
		
		array[0][0] = 4.0;
		array[0][1] = 0.0;
		array[1][0] = 0.0;
		array[1][1] = 2.0;

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		auto evalues = matrix::diagonalize(matrix);

    array = matrix::all_gather(matrix);

		CHECK(real(array[0][0]) == 0.0_a);
		CHECK(imag(array[0][0]) == 0.0_a);
    
		CHECK(real(array[0][1]) == 1.0_a);
		CHECK(imag(array[0][1]) == 0.0_a);
    
		CHECK(real(array[1][0]) == 1.0_a);
		CHECK(imag(array[1][0]) == 0.0_a);
    
		CHECK(real(array[0][0]) == 0.0_a);
		CHECK(imag(array[0][0]) == 0.0_a);
		
    CHECK(evalues[0] == 2.0_a);
		CHECK(evalues[1] == 4.0_a);

	}
	
	SECTION("Real dense 3x3"){
		
		gpu::array<double, 2> array({3, 3});
		
		array[0][0] = 0.088958;
		array[0][1] = 1.183407;
		array[0][2] = 1.191946;
		array[1][0] = 1.183407;
		array[1][1] = 1.371884;
		array[1][2] = 0.705297;
		array[2][0] = 1.191946;
		array[2][1] = 0.705297;
		array[2][2] = 0.392459;

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		auto evalues = matrix::diagonalize(matrix);
		
		CHECK(evalues[0] == -1.0626903983_a);
		CHECK(evalues[1] == 0.1733844724_a);
		CHECK(evalues[2] == 2.7426069258_a);
	}

	SECTION("Complex dense 3x3"){
	
		gpu::array<complex, 2> array({3, 3});
		
		array[0][0] = complex(0.088958,  0.00000);
		array[0][1] = complex(1.183407,  0.08285);
		array[0][2] = complex(1.191946,  0.09413);
		array[1][0] = complex(1.183407, -0.08285);
		array[1][1] = complex(1.371884,  0.00000);
		array[1][2] = complex(0.705297,  0.12840);
		array[2][0] = complex(1.191946, -0.09413);
		array[2][1] = complex(0.705297, -0.12840);
		array[2][2] = complex(0.392459,  0.00000);

		matrix::distributed matrix = matrix::scatter(cart_comm, array, /* root = */ 0);
		
		auto evalues = matrix::diagonalize(matrix);
		
		CHECK(evalues[0] == -1.0703967402_a);
		CHECK(evalues[1] ==  0.1722879629_a);
		CHECK(evalues[2] ==  2.7514097773_a);

	}
}
#endif
