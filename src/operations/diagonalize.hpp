/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__DIAGONALIZE
#define INQ__OPERATIONS__DIAGONALIZE

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa.

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
namespace operations {

template<class Alloc>
auto diagonalize_raw(math::array<double, 2, Alloc>& matrix){

	CALI_CXX_MARK_FUNCTION;

	// the matrix must be square
	assert(std::get<0>(sizes(matrix)) == std::get<1>(sizes(matrix)));

	int nn = std::get<0>(sizes(matrix));
    
	math::array<double, 1> eigenvalues(nn);

	//DATAOPERATIONS RAWLAPACK + CUSOLVER (diagonalization)
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
	dsyev("V", "U", nn, matrix.data_elements(), nn, eigenvalues.data_elements(), &lwork_query, -1, info);

	int lwork = int(lwork_query);
	auto work = (double *) malloc(lwork*sizeof(complex));
    
	dsyev("V", "U", nn, matrix.data_elements(), nn, eigenvalues.data_elements(), work, lwork, info);
	assert(info == 0);
		
	free(work);
#endif

	return eigenvalues;
}

template<class Alloc>
auto diagonalize_raw(math::array<complex, 2, Alloc>& matrix){

	CALI_CXX_MARK_FUNCTION;
	
	// the matrix must be square
	assert(std::get<0>(sizes(matrix)) == std::get<1>(sizes(matrix)));

	int nn = std::get<0>(sizes(matrix));
    
	math::array<double, 1> eigenvalues(nn);

	//DATAOPERATIONS RAWLAPACK + CUSOLVER (diagonalization)
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
	zheev("V", "U", nn, matrix.data_elements(), nn, eigenvalues.data_elements(), &lwork_query, -1, rwork, info);


	int lwork = int(real(lwork_query));
	auto work = (complex *) malloc(lwork*sizeof(complex));
    
	zheev("V", "U", nn, matrix.data_elements(), nn, eigenvalues.data_elements(), work, lwork, rwork, info);
    
	free(rwork);
	free(work);
#endif
		
	return eigenvalues;
}

template <class MatrixType>
auto diagonalize(boost::mpi3::communicator & comm, MatrixType & matrix){

	CALI_CXX_MARK_FUNCTION;
	
	math::array<double, 1> eigenvalues;

	//Dense matrix diagonalization is unstable with respect to small
	//differences in the input. This can cause problems in parallel when
	//each processor runs a copy of the diagonalizer that can give
	//different values. So we diagonalize in one processor and broadcast.
	
	if(comm.rank() == 0){
		eigenvalues = diagonalize_raw(matrix);
	} else {
		eigenvalues = math::array<double, 1>(matrix.size());
	}
	
	if(comm.size() > 1){
		CALI_CXX_MARK_SCOPE("diagonalize::broadcast");

		comm.broadcast_n(raw_pointer_cast(eigenvalues.data_elements()), eigenvalues.num_elements(), 0);
		comm.broadcast_n(raw_pointer_cast(matrix.data_elements()), matrix.num_elements(), 0);
	}

	return eigenvalues;		
}

}
}

///////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_DIAGONALIZE_UNIT_TEST
#undef INQ_OPERATIONS_DIAGONALIZE_UNIT_TEST

#include <catch2/catch.hpp>

#include <operations/randomize.hpp>
#include <math/array.hpp>

TEST_CASE("function operations::diagonalize", "[operations::diagonalize]") {

	auto comm = boost::mpi3::environment::get_world_instance();
	
	SECTION("Real diagonal 2x2"){
	
		using namespace inq;
		using namespace Catch::literals;
		
		math::array<double, 2> matrix({2, 2});
		
		matrix[0][0] = 4.0;
		matrix[0][1] = 0.0;
		matrix[1][0] = 0.0;
		matrix[1][1] = 2.0;
		
		auto evalues = operations::diagonalize(comm, matrix);
		
		CHECK(matrix[0][0] == 0.0_a);
		CHECK(matrix[0][1] == 1.0_a);
		CHECK(matrix[1][0] == 1.0_a);
		CHECK(matrix[0][0] == 0.0_a);
		
		CHECK(evalues[0] == 2.0_a);
		CHECK(evalues[1] == 4.0_a);

	}
	
	SECTION("Complex diagonal 2x2"){
		
		using namespace inq;
		using namespace Catch::literals;
		
		math::array<complex, 2> matrix({2, 2});
		
		matrix[0][0] = 4.0;
		matrix[0][1] = 0.0;
		matrix[1][0] = 0.0;
		matrix[1][1] = 2.0;
		
		auto evalues = operations::diagonalize(comm, matrix);
		
		CHECK(real(matrix[0][0]) == 0.0_a);
		CHECK(imag(matrix[0][0]) == 0.0_a);
		
		CHECK(real(matrix[0][1]) == 1.0_a);
		CHECK(imag(matrix[0][1]) == 0.0_a);
		
		CHECK(real(matrix[1][0]) == 1.0_a);
		CHECK(imag(matrix[1][0]) == 0.0_a);
		
		CHECK(real(matrix[0][0]) == 0.0_a);
		CHECK(imag(matrix[0][0]) == 0.0_a);
		
		CHECK(evalues[0] == 2.0_a);
		CHECK(evalues[1] == 4.0_a);

	}
	
	SECTION("Real dense 3x3"){
	
		using namespace inq;
	using namespace Catch::literals;
		
		math::array<double, 2> matrix({3, 3});
		
		matrix[0][0] = 0.088958;
		matrix[0][1] = 1.183407;
		matrix[0][2] = 1.191946;
		matrix[1][0] = 1.183407;
		matrix[1][1] = 1.371884;
		matrix[1][2] = 0.705297;
		matrix[2][0] = 1.191946;
		matrix[2][1] = 0.705297;
		matrix[2][2] = 0.392459;
		
		auto evalues = operations::diagonalize(comm, matrix);
		
		CHECK(evalues[0] == -1.0626903983_a);
		CHECK(evalues[1] == 0.1733844724_a);
		CHECK(evalues[2] == 2.7426069258_a);
	}

	SECTION("Complex dense 3x3"){
	
		using namespace inq;
		using namespace Catch::literals;
		
		math::array<complex, 2> matrix({3, 3});
		
		matrix[0][0] = complex(0.088958,  0.00000);
		matrix[0][1] = complex(1.183407,  0.08285);
		matrix[0][2] = complex(1.191946,  0.09413);
		matrix[1][0] = complex(1.183407, -0.08285);
		matrix[1][1] = complex(1.371884,  0.00000);
		matrix[1][2] = complex(0.705297,  0.12840);
		matrix[2][0] = complex(1.191946, -0.09413);
		matrix[2][1] = complex(0.705297, -0.12840);
		matrix[2][2] = complex(0.392459,  0.00000);
		
		auto evalues = operations::diagonalize(comm, matrix);
		
		CHECK(evalues[0] == -1.0703967402_a);
		CHECK(evalues[1] ==  0.1722879629_a);
		CHECK(evalues[2] ==  2.7514097773_a);

	}
}



#endif

#endif
