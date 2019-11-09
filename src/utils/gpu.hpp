/* -*- indent-tabs-mode: t -*- */

#ifndef UTILS__GPU
#define UTILS__GPU

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_CUDA
#include <cuda.h>
#endif

#define CUDA_BLOCK_SIZE 1024

#define CUDA_MAX_DIM1 2147483647ULL
#define CUDA_MAX_DIM23 65535

namespace gpu {

	//finds fact1, fact2 < thres such that fact1*fact2 >= val
	void factorize(const size_t val, const size_t thres, size_t & fact1, size_t & fact2){
		fact1 = val;
		fact2 = 1;
		while (fact1 > thres){
			fact1 = (fact1 + 1)/2;
			fact2 *= 2;
		}

		assert(fact1*fact2 >= val);
	}
	
#ifdef HAVE_CUDA
  template <class kernel_type>
  __global__ void cuda_run_kernel_1(unsigned size, kernel_type kernel){
    auto ii = blockIdx.x*blockDim.x + threadIdx.x;
    if(ii < size) kernel(ii);
  }
#endif

  template <class kernel_type>
  void run(long size, kernel_type kernel){

#ifdef HAVE_CUDA

		assert(size <= CUDA_MAX_DIM1);
		
		unsigned nblock = (size + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
    
    cuda_run_kernel_1<<<nblock, CUDA_BLOCK_SIZE>>>(size, kernel);
    
    cudaDeviceSynchronize();
#endif
    
  }

#ifdef HAVE_CUDA
  template <class kernel_type>
  __global__ void cuda_run_kernel_2(unsigned sizex, unsigned sizey, unsigned dim2, kernel_type kernel){
		auto i1 = blockIdx.x*blockDim.x + threadIdx.x;
		auto i2 = blockIdx.y*blockDim.y + threadIdx.y;
		auto i3 = blockIdx.z*blockDim.z + threadIdx.z;
		
    auto ix = i1;
    auto iy = i2 + dim2*i3;
    if(ix < sizex && iy < sizey) kernel(ix, iy);
  }
#endif
 
  template <class kernel_type>
  void run(long sizex, long sizey, kernel_type kernel){

#ifdef HAVE_CUDA
    //OPTIMIZATION, this is not ideal if sizex < CUDA_BLOCK_SIZE
    unsigned nblock = (sizex + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;

		size_t dim2, dim3;
		factorize(sizey, CUDA_MAX_DIM23, dim2, dim3);
		
    cuda_run_kernel_2<<<{nblock, unsigned(dim2), unsigned(dim3)}, {CUDA_BLOCK_SIZE, 1}>>>(sizex, sizey, dim2, kernel);
    
    cudaDeviceSynchronize();
#endif
    
  }

}

#ifdef HAVE_CUDA
#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <math/array.hpp>

namespace gpu {
	
	size_t check_run(size_t size){
		
		math::array<size_t, 1> list(size, 0);

		gpu::run(size,
						 [ll = begin(list)] __device__ (auto ii){
							 atomicAdd((unsigned long long int*) &(ll[ii]), (unsigned long long int) ii + 1);
						 });

		size_t diff = 0;
		for(size_t ii = 0; ii < size; ii++) {
			diff += ii + 1 - list[ii];
		}
		return diff;
	}

	size_t check_run(size_t size1, size_t size2){

		math::array<size_t, 3> list({size1, size2, 2}, 0);

		gpu::run(size1, size2, 
						 [ll = begin(list)] __device__ (auto ii, auto jj){
							 atomicAdd((unsigned long long int*) &(ll[ii][jj][0]), (unsigned long long int) ii + 1);
							 atomicAdd((unsigned long long int*) &(ll[ii][jj][1]), (unsigned long long int) jj + 1);
						 });

		size_t diff = 0;
		for(size_t ii = 0; ii < size1; ii++) {
			for(size_t jj = 0; jj < size2; jj++) {
				diff += ii + 1 - list[ii][jj][0];
				diff += jj + 1 - list[ii][jj][1];
			}
		}

		return diff;

	}
	
}

TEST_CASE("function gpu::run", "[gpu::run]") {

	using namespace Catch::literals;

	SECTION("1D"){
		REQUIRE(gpu::check_run(200) == 0);
		REQUIRE(gpu::check_run(1024) == 0);
		REQUIRE(gpu::check_run(6666) == 0);
		REQUIRE(gpu::check_run(127939) == 0);
	}

	SECTION("2D"){
		REQUIRE(gpu::check_run(200, 200) == 0);
		REQUIRE(gpu::check_run(256, 1200) == 0);
		REQUIRE(gpu::check_run(2023, 4) == 0);
		REQUIRE(gpu::check_run(7, 57*57*57) == 0);		
	}
	
}

#endif
#endif
#endif
