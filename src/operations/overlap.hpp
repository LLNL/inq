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

#include <multi/array.hpp>
#include <basis/field_set.hpp>
#include <cassert>
#include <multi/adaptors/blas.hpp>

#ifdef HAVE_CUDA
#include <multi/memory/adaptors/cuda/allocator.hpp>
#endif

namespace operations {

	template <class field_set_type>
  auto overlap(const field_set_type & phi1, const field_set_type & phi2){

		boost::multi::array<typename field_set_type::value_type, 2>  overlap_matrix({phi1.set_size(), phi1.set_size()});

		boost::multi::blas::gemm('N', 'C', phi1.basis().volume_element(), phi1, phi2, 0.0, overlap_matrix);

		return overlap_matrix;		
  }

	template <class field_set_type>
	auto overlap(const field_set_type & phi){
		return overlap(phi, phi);
	}

#ifdef HAVE_CUDA
	template <class type>
	__global__ void overlap_diagonal_kernel(const long npoints, const int nst, const double vol_element,
																					const type * phi1, const type * phi2, type * overlap){
		
		int ist = blockIdx.x*blockDim.x + threadIdx.x;

		type aa = 0.0;
		for(int ip = 0; ip < npoints; ip++){
			complex p1 = phi1[ip*nst + ist];
			complex p2 = phi2[ip*nst + ist];
			aa += conj(p1)*p2;
		}
		
		overlap[ist] = vol_element*aa;

	}
#endif
	
	template <class field_set_type>
  auto overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

		namespace multi = boost::multi;

		using value_type = typename field_set_type::value_type;
		
		multi::array<value_type, 1>  overlap_vector(phi1.set_size());

		assert(size(overlap_vector) == phi1.set_size());

#ifndef HAVE_CUDA

		//DATAOPERATIONS
		
		//OPTIMIZATION: this can be done more efficiently
    for(int ii = 0; ii < phi1.set_size(); ii++){
			typename field_set_type::value_type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().num_points(); ip++) aa += conj(phi1[ip][ii])*phi2[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
    }

#else

		namespace cuda = multi::memory::cuda;
		
		multi::array<value_type, 1, cuda::allocator<complex>> overlap_cuda(phi1.set_size());
		multi::array<value_type, 2, cuda::allocator<complex>> phi1_cuda = phi1;
		multi::array<value_type, 2, cuda::allocator<complex>> phi2_cuda = phi2;

		//OPTIMIZATION: here we should parallelize over points as well 
		overlap_diagonal_kernel<complex><<<1, phi1.set_size()>>>(phi1.basis().num_points(), phi1.set_size(), phi1.basis().volume_element(),
																														 static_cast<value_type const *>(phi1_cuda.data()),
																														 static_cast<value_type const *>(phi2_cuda.data()),
																														 static_cast<value_type *>(overlap_cuda.data()));

		overlap_vector = overlap_cuda;
		
#endif
		
		return overlap_vector;		
  }
	
	template <class field_set_type>
	auto overlap_diagonal(const field_set_type & phi){
		
		return overlap_diagonal(phi, phi);
	}
	
	template <class field_type>
	auto overlap_single(const field_type & phi1, const field_type & phi2){

		//DATAOPERATIONS
		//OPTIMIZATION: this can be done more efficiently
		typename field_type::value_type overlap = 0.0;
		for(int ip = 0; ip < phi1.basis().num_points(); ip++) overlap += conj(phi1[ip])*phi2[ip];
		return overlap*phi1.basis().volume_element();
	}
	
	template <class field_type>
	auto overlap_single(field_type & phi){
		return overlap_single(phi, phi);
	}
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("function operations::overlap", "[overlap]") {

	using namespace Catch::literals;

}


#endif
#endif
