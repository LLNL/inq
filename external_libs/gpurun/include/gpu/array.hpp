/* -*- indent-tabs-mode: t -*- */

#ifndef GPURUN__MATH__ARRAY
#define GPURUN__MATH__ARRAY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <utils/profiling.hpp>

#include <multi/array.hpp>
// #include <multi/detail/fix_complex_traits.hpp>

#include <complex>

#include <gpu/host.hpp>

#ifdef __NVCC__
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::std::complex<double>> = true;
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::std::complex<float>> = true;
#else  // vvv nvcc (12.1?) doesn't support this kind of customization: "error: expected initializer before ‘<’"
template<class T>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::std::complex<T>> = std::is_trivially_default_constructible<T>::value;
#endif

#ifdef ENABLE_GPU
#include<thrust/complex.h>

#ifdef __NVCC__
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::thrust::complex<double>> = true;
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::thrust::complex<float>> = true;
#else  // vvv nvcc (12.1?) doesn't support this kind of customization: "error: expected initializer before ‘<’"
template<class T>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::thrust::complex<T>> = std::is_trivially_default_constructible<T>::value;
#endif

#ifdef ENABLE_HIP
#define MULTI_USE_HIP
#endif

#include <multi/adaptors/thrust.hpp>
#endif

#ifdef ENABLE_GPU
#include <thrust/mr/disjoint_tls_pool.h>  // for thrust::mr::tls_disjoint_pool
#endif

#ifdef ENABLE_GPU
#include "multi/adaptors/cuda/cublas.hpp" // must be included before blas.hpp
#include "multi/adaptors/cuda/cublas/context.hpp" // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>

namespace gpu {

#ifdef ENABLE_GPU
template<typename Upstream, typename Bookkeeper>
thrust::mr::disjoint_unsynchronized_pool_resource<Upstream, Bookkeeper>& 
LEAKY_tls_disjoint_pool(
			Upstream * upstream = NULL,
			Bookkeeper * bookkeeper = NULL
) {
  static thread_local auto adaptor = new thrust::mr::disjoint_unsynchronized_pool_resource<Upstream, Bookkeeper>(upstream, bookkeeper);
  return *adaptor;
}
template<class T, class Base_ = thrust::mr::allocator<T, thrust::mr::memory_resource<thrust::universal_allocator<void>::pointer>>>
struct caching_allocator : Base_ {
	caching_allocator() : Base_{
		&LEAKY_tls_disjoint_pool(thrust::mr::get_global_resource<thrust::universal_memory_resource>(), thrust::mr::get_global_resource<thrust::mr::new_delete_resource>())
	} {}
	caching_allocator(caching_allocator const&) : caching_allocator{} {}
  	template<class U> struct rebind {using other = caching_allocator<U>;};

  // using Base_::allocate;
  [[nodiscard]] constexpr auto allocate(typename std::allocator_traits<Base_>::size_type n) -> typename std::allocator_traits<Base_>::pointer {
		CALI_MARK_BEGIN("allocate");
    auto ret = std::allocator_traits<Base_>::allocate(*this, n);
		prefetch_to_device(ret, n*sizeof(T), get_current_device());
		CALI_MARK_END("allocate");
    return ret;
  }

  [[nodiscard]] constexpr auto allocate(typename std::allocator_traits<Base_>::size_type n, typename std::allocator_traits<Base_>::const_void_pointer hint) -> typename std::allocator_traits<Base_>::pointer {
		CALI_MARK_BEGIN("allocate_with_hint");
    auto ret = std::allocator_traits<Base_>::allocate(*this, n);
    if(not hint) {
			prefetch_to_device(ret, n*sizeof(T), get_current_device());
      return ret;
    }
		prefetch_to_device(ret, n*sizeof(T), get_device(hint));
		CALI_MARK_END("allocate_with_hint");		
    return ret;
  }
	
	constexpr void deallocate(typename std::allocator_traits<Base_>::pointer p, typename std::allocator_traits<Base_>::size_type n) {
		CALI_MARK_BEGIN("deallocate");
		Base_::deallocate(p, n);	
		CALI_MARK_END("deallocate");	
	}
};
#endif

template <class type, size_t dim,
#ifdef ENABLE_GPU
					class allocator = caching_allocator<type>
#else
					class allocator = std::allocator<type>
#endif
					>
using array = boost::multi::array<type, dim, allocator>;

template <typename ArrayType>
void prefetch(ArrayType const & array){
	prefetch_to_device(array.data_elements(), array.num_elements()*sizeof(typename ArrayType::element_type), get_current_device());
}

template <typename ArrayType>
void prefetch_cpu(ArrayType const & array){
	prefetch_to_device(array.data_elements(), array.num_elements()*sizeof(typename ArrayType::element_type), cpu_device());	
}

}
#endif

#ifdef GPURUN__ARRAY__UNIT_TEST
#undef GPURUN__ARRAY__UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(GPURUN_TEST_FILE, GPURUN_TEST_TAG) {
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
