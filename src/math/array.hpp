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
template<typename Upstream, typename Bookkeeper>
thrust::mr::disjoint_unsynchronized_pool_resource<Upstream, Bookkeeper>& 
LEAKY_tls_disjoint_pool(
			Upstream * upstream = NULL,
			Bookkeeper * bookkeeper = NULL
) {
  static thread_local auto adaptor = new thrust::mr::disjoint_unsynchronized_pool_resource<Upstream, Bookkeeper>(upstream, bookkeeper);
  return *adaptor;
}

template<class T, class Base_ = thrust::mr::allocator<T, thrust::mr::memory_resource<thrust::cuda::universal_pointer<void>>>>
struct caching_allocator : Base_ {
	caching_allocator() : Base_{
		&LEAKY_tls_disjoint_pool(thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(), thrust::mr::get_global_resource<thrust::mr::new_delete_resource>())
	} {}
	caching_allocator(caching_allocator const&) : caching_allocator{} {}
  	template<class U> struct rebind {using other = caching_allocator<U>;};

  // using Base_::allocate;
  [[nodiscard]] constexpr auto allocate(typename std::allocator_traits<Base_>::size_type n) -> typename std::allocator_traits<Base_>::pointer {
		CALI_CXX_MARK_SCOPE("allocate");
    auto ret = std::allocator_traits<Base_>::allocate(*this, n);
		prefetch_to_device(ret, n*sizeof(T), get_current_device());
    return ret;
  }

  [[nodiscard]] constexpr auto allocate(typename std::allocator_traits<Base_>::size_type n, typename std::allocator_traits<Base_>::const_void_pointer hint) -> typename std::allocator_traits<Base_>::pointer {
		CALI_CXX_MARK_SCOPE("allocate_with_hint");
    auto ret = std::allocator_traits<Base_>::allocate(*this, n);
    if(not hint) {
			prefetch_to_device(ret, n*sizeof(T), get_current_device());
      return ret;
    }
		prefetch_to_device(ret, n*sizeof(T), get_device(hint));
    return ret;
  }
	
	constexpr void deallocate(typename std::allocator_traits<Base_>::pointer p, typename std::allocator_traits<Base_>::size_type n) {
		CALI_CXX_MARK_SCOPE("deallocate");
		Base_::deallocate(p, n);
	}
	
private:
  using device_index = int;
  static auto get_current_device() -> device_index {
    int device;
    switch(cudaGetDevice(&device)) {
    case cudaSuccess          : break;
    case cudaErrorInvalidValue: assert(0);
    }
    return device;
  }
  static void prefetch_to_device(typename std::allocator_traits<Base_>::const_void_pointer p, typename std::allocator_traits<Base_>::size_type byte_count, device_index d) {
    switch(cudaMemPrefetchAsync(raw_pointer_cast(p), byte_count, d)) {
    case cudaSuccess           : return;
    case cudaErrorInvalidValue : ;
    case cudaErrorInvalidDevice: ;
    }
    assert(0);
  }

  static auto get_device(typename std::allocator_traits<Base_>::const_void_pointer p) -> device_index {
    cudaPointerAttributes attr{};
    switch(cudaPointerGetAttributes(&attr, raw_pointer_cast(p))) {
    case cudaSuccess           : break;
    case cudaErrorInvalidDevice:
    case cudaErrorInvalidValue : assert(0);
    }
    assert(attr.type == cudaMemoryTypeManaged);
    return attr.device;
  }
};
#endif

template <class type, size_t dim,
#ifdef ENABLE_CUDA
					class allocator = caching_allocator<type>
#else
					class allocator = std::allocator<type>
#endif
					>
using array = boost::multi::array<type, dim, allocator>;

template <class type, size_t dim,
#ifdef ENABLE_CUDA
					class allocator = caching_allocator<type>
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
	int device;
	cudaGetDevice(&device);
	cudaMemPrefetchAsync(raw_pointer_cast(array.data_elements()), array.num_elements()*sizeof(typename ArrayType::element_type), device);
#endif
}

template <typename ArrayType>
void prefetch_cpu(ArrayType const &
#ifdef ENABLE_CUDA
							array
#endif
							){
#ifdef ENABLE_CUDA
	cudaMemPrefetchAsync(raw_pointer_cast(array.data_elements()), array.num_elements()*sizeof(typename ArrayType::element_type), cudaCpuDeviceId);
#endif
}

}
}
#endif

#ifdef INQ_MATH_ARRAY_UNIT_TEST
#undef INQ_MATH_ARRAY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
