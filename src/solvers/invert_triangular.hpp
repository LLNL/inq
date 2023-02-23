/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__INVERT_TRIANGULAR
#define INQ__SOLVERS__INVERT_TRIANGULAR

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

#include <FC.h>

#include <tuple> //std::get
#include <cassert>

#include <math/array.hpp>
#include <math/complex.hpp>
#include <math/subspace_matrix.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <cusolverDn.h>

#define dtrtri FC_GLOBAL(dtrtri, DTRTRI) 
extern "C" void  dtrtri(const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);

#define ztrtri FC_GLOBAL(ztrtri, ZTRTRI) 
extern "C" void  ztrtri(const char * uplo, const char * diag, const int * n, inq::complex * a, const int * lda, int * info);

namespace inq {
namespace solvers {

template <typename Type>
void invert_triangular(math::subspace_matrix<Type> & matrix){
  CALI_CXX_MARK_SCOPE("invert_triangular_double");

	static_assert(std::is_same_v<Type, double> or std::is_same_v<Type, complex>, "invert_triangular is only implemented for double and complex");
	
	int nn = std::get<0>(sizes(matrix.array()));

	auto matrix_data = raw_pointer_cast(matrix.array().data_elements());
	int info;
	
#ifndef ENABLE_CUDA	

	if constexpr (std::is_same_v<Type, double>){
		dtrtri("U", "N", &nn, matrix_data, &nn, &info);
	} else {
		ztrtri("U", "N", &nn, matrix_data, &nn, &info);
	}

#else

	cudaDataType datatype = CUDA_C_64F;
	if constexpr (std::is_same_v<Type, double>){
		datatype = CUDA_R_64F;
	}
	
	std::cout << CUSOLVER_STATUS_SUCCESS << '\t' << CUSOLVER_STATUS_NOT_INITIALIZED << '\t' << CUSOLVER_STATUS_NOT_SUPPORTED << '\t' << CUSOLVER_STATUS_INVALID_VALUE << '\t' << CUSOLVER_STATUS_INTERNAL_ERROR << std::endl;
	
	cusolverDnHandle_t cusolver_handle;
	
	[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	
	size_t worksize_host, worksize_device;

	cusolver_status = cusolverDnXtrtri_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, nn, datatype, matrix_data, nn, &worksize_device, &worksize_host);
	std::cout << cusolver_status << std::endl;
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

	std::cout << worksize_device << '\t' << worksize_host << std::endl;

	void * work_device;
	cudaMalloc(&work_device, worksize_device);
	//math::array<char, 1> work_device(worksize_device);
	//	assert(work_device.num_elements() == worksize_device);

	auto work_host = malloc(worksize_host);
	assert(work_host != NULL);
	
	cusolver_status = cusolverDnXtrtri(cusolver_handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, nn, datatype, matrix_data, nn, work_device, worksize_device, work_host, worksize_host, &info);
	std::cout << cusolver_status << '\t' << info << std::endl;	
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

	free(work_host);
	cusolverDnDestroy(cusolver_handle);
	
#endif
	
	gpu::run(nn, nn, [mat = begin(matrix.array())] GPU_LAMBDA (auto ii, auto jj){
		if(ii < jj) mat[ii][jj] = 0.0;
	});
}

}
}
#endif

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_INVERT_TRIANGULAR_UNIT_TEST
#undef INQ_SOLVERS_INVERT_TRIANGULAR_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <mpi3/environment.hpp>

#include <math/array.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  auto comm = boost::mpi3::environment::get_world_instance();
	inq::parallel::cartesian_communicator<2> cart_comm(comm, {});

	SECTION("Real 2x2"){
	
		using namespace inq;
		using namespace Catch::literals;

    math::subspace_matrix<double> matrix(cart_comm, 2);
		
		matrix.array()[0][0] = 4.0;
		matrix.array()[1][0] = -1.0;
		matrix.array()[1][1] = 2.0;

		solvers::invert_triangular(matrix);
    
    CHECK(matrix.array()[0][0] == 0.25);
		CHECK(matrix.array()[0][1] == 0.0);
    CHECK(matrix.array()[1][0] == 0.125);
    CHECK(matrix.array()[1][1] == 0.5);
  }

  SECTION("Complex 2x2"){
	
		using namespace inq;
		using namespace Catch::literals;

    math::subspace_matrix<complex> matrix(cart_comm, 2);
		
		matrix.array()[0][0] = 4.0;
		matrix.array()[1][0] = -1.0;
		matrix.array()[1][1] = 2.0;

		solvers::invert_triangular(matrix);
    
    CHECK(matrix.array()[0][0] == 0.25);
		CHECK(matrix.array()[0][1] == 0.0);
    CHECK(matrix.array()[1][0] == 0.125);
    CHECK(matrix.array()[1][1] == 0.5);
  }
	
  SECTION("Complex 3x3"){
	
		using namespace inq;
		using namespace Catch::literals;

    math::subspace_matrix<complex> matrix(cart_comm, 3);
		
		matrix.array()[0][0] = 4.0;
		matrix.array()[0][1] = 0.0;
		matrix.array()[0][2] = 0.0;
		matrix.array()[1][0] = 0.0;
		matrix.array()[1][1] = 10.0;
		matrix.array()[1][2] = 0.0;
		matrix.array()[2][0] = 0.0;
		matrix.array()[2][1] = 0.0;
		matrix.array()[2][2] = 3.0;		

		solvers::invert_triangular(matrix);

  }
	
	SECTION("Complex NxN"){
	
		using namespace inq;
		using namespace Catch::literals;

		auto nn = 100;

    math::subspace_matrix<complex> matrix(cart_comm, nn);
		
		for(int ii = 0; ii < nn; ii++){
			for(int jj = 0; jj < nn; jj++){
				matrix.array()[ii][jj] = (ii == jj) ? 1.0 : 0.0;
			}
		}
		
		solvers::invert_triangular(matrix);
		
  }
}
#endif
