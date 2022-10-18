/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__ARRAY
#define INQ__MATH__ARRAY

/*
 Copyright (C) 2019-2022 Xavier Andrade, Alfredo A. Correa.

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

#include <utils/profiling.hpp>

#include <multi/detail/fix_complex_traits.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/thrust/fix_complex_traits.hpp>
#include <multi/adaptors/thrust.hpp>
#endif

#include <multi/array.hpp>

#ifdef ENABLE_CUDA
#include <thrust/system/cuda/memory.h>  // for ::thrust::cuda::universal_allocator<type>
#include <thrust/mr/disjoint_tls_pool.h>  // for thrust::mr::tls_disjoint_pool
#endif

#ifdef ENABLE_CUDA
#include "multi/adaptors/cuda/cublas.hpp" // must be included before blas.hpp
#include "multi/adaptors/cuda/cublas/context.hpp" // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>

namespace inq {
namespace math {

#ifdef ENABLE_CUDA
template<class T, class Base_ = thrust::mr::allocator<T, thrust::mr::memory_resource<thrust::cuda::universal_pointer<void>>>>
struct caching_allocator : Base_ {
	caching_allocator() : Base_{
		&thrust::mr::tls_disjoint_pool(thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(), thrust::mr::get_global_resource<thrust::mr::new_delete_resource>())
	} {}
	caching_allocator(caching_allocator const&) : caching_allocator{} {}
	template<class U> struct rebind {using other = caching_allocator<U>;};
};
#endif

template <class type, size_t dim,
#ifdef ENABLE_CUDA
					class allocator = caching_allocator<type>  // ::thrust::cuda::universal_allocator<type>
#else
					class allocator = std::allocator<type>
#endif
					>
using array = boost::multi::array<type, dim, allocator>;

template <class type, size_t dim,
#ifdef ENABLE_CUDA
					class allocator = ::thrust::cuda::universal_allocator<type>
#else
					class allocator = std::allocator<type>
#endif
					>
using array_nopre = boost::multi::array<type, dim, allocator>;

template <typename ArrayType>
void prefetch(ArrayType const &
#ifdef ENABLE_CUDA
							array
#endif
							){
#ifdef ENABLE_CUDA
	cudaMemPrefetchAsync(raw_pointer_cast(array.data_elements()), array.num_elements()*sizeof(typename ArrayType::element_type), 0);
#endif
}


}
}

#ifdef INQ_MATH_ARRAY_UNIT_TEST
#undef INQ_MATH_ARRAY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("math::array", "[math::array]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
