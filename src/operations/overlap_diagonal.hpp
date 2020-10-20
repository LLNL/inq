/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP_DIAGONAL
#define OPERATIONS__OVERLAP_DIAGONAL

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

#include <inq_config.h>

#include <math/array.hpp>
#include <cassert>
#ifdef ENABLE_CUDA
#include "multi/adaptors/blas/cuda.hpp" // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>
#include <operations/integral.hpp>

#include <caliper/cali.h>

#include <gpu/run.hpp>

namespace inq {
namespace operations {

template <class phi1p_type, class phi2p_type, class overlap_type>
__global__ void overlap_diagonal_kernel(const long npoints, const long nst, const double vol_element, const phi1p_type phi1p, const phi2p_type phi2p, overlap_type overlap){

	unsigned int ist = blockIdx.x*blockDim.x + threadIdx.x;

	if(ist >= nst) return;
	
	typename overlap_type::element aa = 0.0;
	
	for(int ip = 0; ip < npoints; ip++){
		auto p1 = phi1p[ip][ist];
		auto p2 = phi2p[ip][ist];
		aa += conj(p1)*p2;
	}

	overlap[ist] = vol_element*aa;

}

template <class field_set_type>
math::array<typename field_set_type::element_type, 1> overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

	CALI_CXX_MARK_SCOPE("overlap_diagonal 2 arg");
	
	using type = typename field_set_type::element_type;
		
	math::array<type, 1> overlap_vector(phi1.set_part().local_size());

	assert(size(overlap_vector) == phi2.set_part().local_size());

	if(phi2.set_part().local_size() == 1){
#ifdef ENABLE_CUDA
		if(typeid(typename field_set_type::element_type) == typeid(complex)) {
			cublasZdotc(boost::multi::cublas::global_context().get(), phi1.basis().part().local_size(),
									(const cuDoubleComplex *) raw_pointer_cast(phi1.matrix().data_elements()), 1, (const cuDoubleComplex *)  raw_pointer_cast(phi2.matrix().data_elements()), 1,
									(cuDoubleComplex *) raw_pointer_cast(overlap_vector.data_elements()));
		} else {
			cublasDdot(boost::multi::cublas::global_context().get(), phi1.basis().part().local_size(),
								 (const double *) raw_pointer_cast(phi1.matrix().data_elements()), 1, (const double *) raw_pointer_cast(phi2.matrix().data_elements()), 1,
								 (double *) raw_pointer_cast(overlap_vector.data_elements()));
		}
#else

		using boost::multi::blas::dot;
		using boost::multi::blas::conj;
		
		overlap_vector[0] = dot(boost::multi::blas::conj(phi1.matrix().rotated()[0]), phi2.matrix().rotated()[0]);
#endif
		overlap_vector[0] *= phi1.basis().volume_element();
	} else {
		
		//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifndef ENABLE_CUDA
		
		//OPTIMIZATION: this can be done more efficiently
		for(int ii = 0; ii < phi1.local_set_size(); ii++){
			type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().part().local_size(); ip++) aa += conj(phi1.matrix()[ip][ii])*phi2.matrix()[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
		}
		
#else

		int mingridsize = 0;
		int blocksize = 0;
		gpu::check_error(cudaOccupancyMaxPotentialBlockSize(&mingridsize, &blocksize,  overlap_diagonal_kernel<decltype(begin(phi1.matrix())), decltype(begin(phi2.matrix())), decltype(begin(overlap_vector))>));
	
		auto size = phi1.local_set_size();
		unsigned nblock = (size + blocksize - 1)/blocksize;
		
		overlap_diagonal_kernel<<<nblock, blocksize>>>(phi1.basis().part().local_size(), phi1.local_set_size(), phi1.basis().volume_element(), begin(phi1.matrix()), begin(phi2.matrix()), begin(overlap_vector));

		gpu::check_error(cudaGetLastError());    
		
		cudaDeviceSynchronize();
		
#endif

	}
	
	if(phi1.basis().part().parallel()){
		phi1.basis().comm().all_reduce_in_place_n(static_cast<type *>(overlap_vector.data()), overlap_vector.size(), std::plus<>{});
	}
		
	return overlap_vector;
}
	
template <class field_set_type>
auto overlap_diagonal(const field_set_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_diagonal 1 arg");
	
	return overlap_diagonal(phi, phi);
}
	
}
}


#ifdef INQ_OPERATIONS_OVERLAP_DIAGONAL_UNIT_TEST
#undef INQ_OPERATIONS_OVERLAP_DIAGONAL_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::overlap_diagonal", "[operations::overlap_diagonal]") {
	
	using namespace inq;
	using namespace Catch::literals;

	const int npoint = 100;
	const int nvec = 12;
			
	auto comm = boost::mpi3::environment::get_world_instance();
		
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {1, comm.size()});

	auto basis_comm = cart_comm.axis(1);

	CHECK(basis_comm.size() == comm.size());
		
	basis::trivial bas(npoint, basis_comm);

	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, double> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig + 1)*sqrt(jjg);
				bb.matrix()[ii][jj] = -0.05/(iig + 1)*sqrt(jjg);
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa, bb);
				
			for(int jj = 0; jj < nvec; jj++) CHECK(dd[jj] == Approx(-jj));
		}
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig)*sqrt(jjg);
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa);
								
			for(int jj = 0; jj < nvec; jj++) CHECK(dd[jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
		}
					
			
	}

	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig + 1)*sqrt(jjg)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig));
				bb.matrix()[ii][jj] = -0.05/(iig + 1)*sqrt(jjg)*exp(complex(0.0, M_PI/4 + M_PI/7*iig));
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa, bb);

			CHECK(std::get<0>(sizes(dd)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++){
				CHECK(fabs(real(dd[jj])) < 1.0e-14);
				CHECK(imag(dd[jj]) == Approx(-jj));
			}
		}
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig)*sqrt(jjg)*exp(complex(0.0, M_PI/65.0*iig));
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa);

			CHECK(std::get<0>(sizes(dd)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++) CHECK(real(dd[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
		}
					
	}

}

#endif
#endif
