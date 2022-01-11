/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__ORTHOGONALIZE
#define INQ__OPERATIONS__ORTHOGONALIZE

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa

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

#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>
#include <operations/overlap.hpp>

#ifdef ENABLE_CUDA
#include <cusolverDn.h>
#endif

#include "FC.h"

#define zpotrf FC_GLOBAL(zpotrf, ZPOTRF) 
extern "C" void zpotrf(const char * uplo, const int * n, inq::complex * a, const int * lda, int * info);

//#define blas_ztrsm FC_GLOBAL(ztrsm, ZTRSM) 
//extern "C" void blas_ztrsm(const char& side, const char& uplo, const char& transa, const char& diag,
//											const long& m, const long& n, const complex& alpha, const complex * a, const long& lda, complex * B, const long& ldb);

namespace inq {
namespace operations {

template <class field_set_type>
void orthogonalize(field_set_type & phi, bool nocheck = false){

	CALI_CXX_MARK_FUNCTION;

	assert(phi.set_comm().size() == 1);
	
	auto olap = overlap(phi);

	const int nst = phi.set_size();
	int info;
	
	//DATAOPERATIONS RAWLAPACK zpotrf
#ifdef ENABLE_CUDA
	{

		CALI_CXX_MARK_SCOPE("cuda_zpotrf");
		
		cusolverDnHandle_t cusolver_handle;
			
		[[maybe_unused]] auto cusolver_status = cusolverDnCreate(&cusolver_handle);
		assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
			
		//query the work size
		int lwork;
		cusolver_status = cusolverDnZpotrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(olap.data_elements()), nst, &lwork);
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

		cusolver_status = cusolverDnZpotrf(cusolver_handle, CUBLAS_FILL_MODE_UPPER, nst, (cuDoubleComplex *) raw_pointer_cast(olap.data_elements()), nst, work, lwork, devInfo);
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
		zpotrf("U", &nst, raw_pointer_cast(olap.data_elements()), &nst, &info);
	}
#endif

	if(not nocheck and info < 0){
		std::printf("Error: Failed orthogonalization in ZPOTRF! info is %10i.\n", info);
		abort();
	} else if(info > 0) {
		std::printf("Warning: Imperfect orthogonalization in ZPOTRF! info is %10i, subspace size is %10i\n", info, nst); 
	}

	{
		CALI_CXX_MARK_SCOPE("orthogonalize_trsm");
		namespace blas = boost::multi::blas;
		blas::trsm(blas::side::right, blas::filling::upper, 1., blas::H(olap), phi.matrix());
	}
}

template <class FieldSetType1, class FieldSetType2>
void orthogonalize_single(FieldSetType1 & vec, FieldSetType2 const & phi, int num_states = -1){

	if(num_states == 0) return;
	
	CALI_CXX_MARK_FUNCTION;

	if(num_states == -1) num_states = phi.set_size();
		
	assert(num_states <= phi.set_size());

	namespace blas = boost::multi::blas;

	//the matrix of phi restricted to the vectors we are going to use
	auto phi_restricted = phi.matrix()({0, phi.basis().local_size()}, {0, num_states});

	math::array<typename FieldSetType1::element_type, 2> olap;

	if(num_states == 1){
		// avoid a bug in multi by making an explicit copy
		//   https://gitlab.com/correaa/boost-multi/-/issues/97
		math::array<typename FieldSetType1::element_type, 2> phi_restricted_copy = phi_restricted;
		olap = blas::gemm(phi.basis().volume_element(), blas::H(phi_restricted_copy), vec.matrix());
	} else {
		//this should be done by gemv, but while multi gets better gemv support we use gemm
		olap = blas::gemm(phi.basis().volume_element(), blas::H(phi_restricted), vec.matrix());
	}
	
	if(phi.basis().comm().size() > 1){
		phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(olap.data_elements()), olap.num_elements(), std::plus<>{});
	}

	vec.matrix() += blas::gemm(-1.0, phi_restricted, olap);

}

}
}

#ifdef INQ_OPERATIONS_ORTHOGONALIZE_UNIT_TEST
#undef INQ_OPERATIONS_ORTHOGONALIZE_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::orthogonalize", "[operations::orthogonalize]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using math::vector3;

	systems::box box = systems::box::cubic(6.3_b).cutoff_energy(25.0_Ha);
	basis::real_space pw(box);

	SECTION("Dimension 3"){
		basis::field_set<basis::real_space, complex> phi(pw, 3);
		
		operations::randomize(phi);
		
		operations::orthogonalize(phi);
		
		auto olap = operations::overlap(phi);
		
		std::cout << "------" << std::endl;
		
		std::cout << olap[0][0] << '\t' << olap[0][1] << '\t' << olap[0][2] << std::endl;
		std::cout << olap[1][0] << '\t' << olap[1][1] << '\t' << olap[1][2] << std::endl;
		std::cout << olap[2][0] << '\t' << olap[2][1] << '\t' << olap[2][2] << std::endl;
		
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					CHECK(real(olap[ii][ii]) == 1.0_a);
					CHECK(fabs(imag(olap[ii][ii])) < 1e-14);
			} else {
					CHECK(fabs(olap[ii][jj]) < 1e-14);
				}
			}
		}
	}
	
	SECTION("Dimension 2, warning test"){
		basis::field_set<basis::real_space, complex> phi(pw, 2);

		operations::randomize(phi);
		
		int nbasis = std::get<0>(sizes(phi.matrix()));
		
	        for(int ii = 0 ; ii < phi.set_size(); ii++){
		  for(int jj = 0; jj < nbasis; jj++){
		    phi.matrix().transposed()[ii][jj]=0.0;
		    //if(ii == jj) phi.matrix().transposed()[ii][jj] = 1.0;
		  }
		}
		phi.matrix().transposed()[0][0] = 1.0;
		phi.matrix().transposed()[1][0] = 1.0;
		
		operations::orthogonalize(phi, /* ncheck = */ true);
		
		//auto olap = operations::overlap(phi);
		
		std::cout << "------" << std::endl;
		
		
	}

	SECTION("Dimension 100"){
		basis::field_set<basis::real_space, complex> phi(pw, 100);
		
		operations::randomize(phi);
		
		operations::orthogonalize(phi);
		
		auto olap = operations::overlap(phi);
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					CHECK(real(olap[ii][ii]) == 1.0_a);
					CHECK(fabs(imag(olap[ii][ii])) < 1e-14);
				} else {
					CHECK(fabs(olap[ii][jj]) < 1e-13);
				}
			}
		}
	}


	SECTION("Dimension 37 - double orthogonalize"){
		basis::field_set<basis::real_space, complex> phi(pw, 37);
		
		operations::randomize(phi);
		
		operations::orthogonalize(phi);
		operations::orthogonalize(phi);
		
		auto olap = operations::overlap(phi);
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					CHECK(real(olap[ii][ii]) == 1.0_a);
					CHECK(fabs(imag(olap[ii][ii])) < 1e-16);
				} else {
					CHECK(fabs(olap[ii][jj]) < 5e-16);
				}
			}
		}
	}
	
	SECTION("single -- Dimension 3"){
		basis::field_set<basis::real_space, complex> phi(pw, 3);
		basis::field_set<basis::real_space, complex> vec(pw, 1);
		
		operations::randomize(phi);
		operations::orthogonalize(phi);

		operations::randomize(vec);
		operations::orthogonalize_single(vec, phi, 1);

		for(int ist = 0; ist < 1; ist++){
			
			complex olap = 0.0;
			
			for(long ip = 0; ip < phi.basis().local_size(); ip++){
				olap += conj(phi.matrix()[ip][ist])*vec.matrix()[ip][0];
			}

			phi.basis().comm().all_reduce_in_place_n(&olap, 1, std::plus<>{});
	
			CHECK(fabs(olap) < 5e-14);
		}
		
		operations::randomize(vec);
		operations::orthogonalize_single(vec, phi, 2);

		for(int ist = 0; ist < 2; ist++){
			
			complex olap = 0.0;
			
			for(long ip = 0; ip < phi.basis().local_size(); ip++){
				olap += conj(phi.matrix()[ip][ist])*vec.matrix()[ip][0];
			}

			phi.basis().comm().all_reduce_in_place_n(&olap, 1, std::plus<>{});
	
			CHECK(fabs(olap) < 5e-14);
		}
		
		operations::randomize(vec);
		operations::orthogonalize_single(vec, phi);

		for(int ist = 0; ist < 3; ist++){
			
			complex olap = 0.0;
			
			for(long ip = 0; ip < phi.basis().local_size(); ip++){
				olap += conj(phi.matrix()[ip][ist])*vec.matrix()[ip][0];
			}

			phi.basis().comm().all_reduce_in_place_n(&olap, 1, std::plus<>{});
	
			CHECK(fabs(olap) < 5e-14);
		}
		
	}
}


#endif

#endif
