/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GPU__REDUCE
#define INQ__GPU__REDUCE

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

#include <inq_config.h>

#ifdef ENABLE_CUDA
#include <cuda.h>
#endif

#include <cassert>

#include <gpu/run.hpp>
#include <math/array.hpp>

namespace inq {
namespace gpu {

#ifdef ENABLE_CUDA
template <class kernel_type, class array_type>
__global__ void reduce_kernel_1d(long size, kernel_type kernel, array_type odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename array_type::element *) shared_mem;
	
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;
	unsigned int ii = blockIdx.x*blockDim.x + threadIdx.x;

	if(ii < size){
		reduction_buffer[tid] = kernel(ii);
	} else {
		reduction_buffer[tid] = (typename array_type::element) 0.0;
	}

	__syncthreads();

	// do reduction in shared mem
	for (unsigned int s = blockDim.x/2; s > 0; s >>= 1){
		if (tid < s) {
			reduction_buffer[tid] += reduction_buffer[tid + s];
		}
		__syncthreads();
	}
	
	// write result for this block to global mem
	if (tid == 0) odata[blockIdx.x] = reduction_buffer[0];

}
#endif

template <typename array_type>
struct array_access {
  array_type array;

  GPU_FUNCTION auto operator()(long ii) const {
    return array[ii];
  }

  GPU_FUNCTION auto operator()(long ix, long iy) const {
    return array[ix][iy];
  }
  
};

template <class kernel_type>
auto reduce(long size, kernel_type kernel) -> decltype(kernel(0)) {

  using type = decltype(kernel(0));
  
#ifndef ENABLE_CUDA

  type accumulator = 0.0;
  for(long ii = 0; ii < size; ii++){
    accumulator += kernel(ii);
  }
  return accumulator;

#else

	const int blocksize = 1024;

	unsigned nblock = (size + blocksize - 1)/blocksize;
	math::array<type, 1> result(nblock);

  reduce_kernel_1d<<<nblock, blocksize, blocksize*sizeof(type)>>>(size, kernel, begin(result));	
  check_error(cudaGetLastError());
	
  if(nblock == 1) {
    cudaDeviceSynchronize();
    return result[0];
  } else {
    return reduce(nblock, array_access<decltype(begin(result))>{begin(result)});
  }
  
#endif
}

#ifdef ENABLE_CUDA
template <class kernel_type, class array_type>
__global__ void reduce_kernel_2d(long sizex, long sizey, kernel_type kernel, array_type odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename array_type::element *) shared_mem; // {blockDim.x, blockDim.y}
	
	// each thread loads one element from global to shared mem
  unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.y;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;

	if(iy < sizey){
		reduction_buffer[threadIdx.x + blockDim.x*tid] = kernel(ix, iy);
	} else {
		reduction_buffer[threadIdx.x + blockDim.x*tid] = (typename array_type::element) 0.0;
	}

	__syncthreads();

	// do reduction in shared mem
	for (unsigned int s = blockDim.y/2; s > 0; s >>= 1){
		if (tid < s) {
			reduction_buffer[threadIdx.x + blockDim.x*tid] += reduction_buffer[threadIdx.x + blockDim.x*(tid + s)];
		}
		__syncthreads();
	}
	
	// write result for this block to global mem
	if (tid == 0) odata[blockIdx.y][ix] = reduction_buffer[threadIdx.x];

}
#endif


template <class kernel_type>
auto reduce(long sizex, long sizey, kernel_type kernel) -> math::array<decltype(kernel(0, 0)), 1> {

  using type = decltype(kernel(0, 0));

#ifndef ENABLE_CUDA

  math::array<type, 1> accumulator(sizex, 0.0);

  for(long iy = 0; iy < sizey; iy++){
    for(long ix = 0; ix < sizex; ix++){
      accumulator[ix] += kernel(ix, iy);
    }
  }
  
  return accumulator;
  
#else

  const int blocksize = 1024;

	unsigned nblock = (sizey + blocksize - 1)/blocksize;
	math::array<type, 2> result({nblock, sizex}, 666.0);

	struct dim3 dg{(unsigned) sizex, nblock};
  struct dim3 db{1, blocksize};

  auto shared_mem_size = blocksize*sizeof(type);

  assert(shared_mem_size <= 48*1024);
  
  reduce_kernel_2d<<<dg, db, shared_mem_size>>>(sizex, sizey, kernel, begin(result));	
  check_error(cudaGetLastError());
	
  if(nblock == 1) {
    cudaDeviceSynchronize();

		assert(result[0].size() == sizex);
		
    return result[0];
  } else {
    return reduce(sizex, nblock, array_access<decltype(begin(result.transposed()))>{begin(result.transposed())});		
  }
  
#endif

}

}
}

#ifdef INQ_GPU_REDUCE_UNIT_TEST
#undef INQ_GPU_REDUCE_UNIT_TEST

#include <catch2/catch.hpp>

struct ident {
  GPU_FUNCTION auto operator()(long ii) const {
    return double(ii);
  }
};

struct prod {
  GPU_FUNCTION auto operator()(long ix, long iy) const {
    return double(ix)*double(iy);
  }
};

TEST_CASE("function gpu::reduce", "[gpu::reduce]") {
  
	using namespace inq;
	using namespace Catch::literals;

  SECTION("1D"){
    const long maxsize = 129140163;
    
    for(long nn = 1; nn <= maxsize; nn *= 3){
			CHECK(gpu::reduce(nn, ident{}) == (nn*(nn - 1.0)/2.0));    
    }
  }

  SECTION("2D"){
    
    const long maxsize = 390625;
    
    for(long nx = 1; nx <= 10000; nx *= 10){
      for(long ny = 1; ny <= maxsize; ny *= 5){
        
        auto res = gpu::reduce(nx, ny, prod{});

				CHECK(typeid(decltype(res)) == typeid(math::array<double, 1>));
				
        CHECK(res.size() == nx);
        for(long ix = 0; ix < nx; ix++) CHECK(res[ix] == double(ix)*ny*(ny - 1.0)/2.0);
      }
    }
  }
  
}

#endif
#endif

