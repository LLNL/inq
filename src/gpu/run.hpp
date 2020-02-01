/* -*- indent-tabs-mode: t -*- */

#ifndef GPU__RUN
#define GPU__RUN

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

#ifdef HAVE_CUDA
#define GPU_FUNCTION __host__ __device__
#else
#define GPU_FUNCTION
#endif

#define CUDA_BLOCK_SIZE 1024

#define CUDA_MAX_DIM1 2147483647ULL
#define CUDA_MAX_DIM23 65535

namespace gpu {

	//finds fact1, fact2 < thres such that fact1*fact2 >= val
	inline static void factorize(const size_t val, const size_t thres, size_t & fact1, size_t & fact2){
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
  void run(size_t size, kernel_type kernel){

#ifdef HAVE_CUDA

		assert(size <= CUDA_MAX_DIM1);
		
		unsigned nblock = (size + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
    
    cuda_run_kernel_1<<<nblock, CUDA_BLOCK_SIZE>>>(size, kernel);

		assert(cudaGetLastError() == CUDA_SUCCESS);
		
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
  void run(size_t sizex, size_t sizey, kernel_type kernel){

#ifdef HAVE_CUDA
    //OPTIMIZATION, this is not ideal if sizex < CUDA_BLOCK_SIZE
    unsigned nblock = (sizex + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;

		size_t dim2, dim3;
		factorize(sizey, CUDA_MAX_DIM23, dim2, dim3);

		struct dim3 dg{nblock, unsigned(dim2), unsigned(dim3)};
		struct dim3 db{CUDA_BLOCK_SIZE, 1, 1};
    cuda_run_kernel_2<<<dg, db>>>(sizex, sizey, dim2, kernel);
    
		assert(cudaGetLastError() == CUDA_SUCCESS);
		
    cudaDeviceSynchronize();
#endif
    
  }
	
#ifdef HAVE_CUDA
  template <class kernel_type>
  __global__ void cuda_run_kernel_3(unsigned sizex, unsigned sizey, unsigned sizez, kernel_type kernel){
    auto ix = blockIdx.x*blockDim.x + threadIdx.x;
    auto iy = blockIdx.y*blockDim.y + threadIdx.y;
		auto iz = blockIdx.z*blockDim.z + threadIdx.z;
    if(ix < sizex && iy < sizey && iz < sizez) kernel(ix, iy, iz);

  }
#endif
 
  template <class kernel_type>
  void run(size_t sizex, size_t sizey, size_t sizez, kernel_type kernel){

#ifdef HAVE_CUDA
    //OPTIMIZATION, this is not ideal if sizex < CUDA_BLOCK_SIZE
    unsigned nblock = (sizex + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
    struct dim3 dg{nblock, unsigned(sizey), unsigned(sizez)};
		struct dim3 db{CUDA_BLOCK_SIZE, 1, 1};
		cuda_run_kernel_3<<<dg, db>>>(sizex, sizey, sizez, kernel);

		assert(cudaGetLastError() == CUDA_SUCCESS);
		
    cudaDeviceSynchronize();
#endif
    
  }
	
#ifdef HAVE_CUDA
  template <class kernel_type>
  __global__ void cuda_run_kernel_4(unsigned sizex, unsigned sizey, unsigned sizez, unsigned sizew, kernel_type kernel){
    auto ix = blockIdx.x*blockDim.x + threadIdx.x;
    auto iy = blockIdx.y*blockDim.y + threadIdx.y;
		auto iz = blockIdx.z*blockDim.z + threadIdx.z;
    if(ix < sizex && iy < sizey && iz < sizez){
			for(int iw = 0; iw < sizew; iw++){
				kernel(ix, iy, iz, iw);
			}
		}
  }
#endif
 
  template <class kernel_type>
  void run(size_t sizex, size_t sizey, size_t sizez, size_t sizew, kernel_type kernel){

#ifdef HAVE_CUDA
    //OPTIMIZATION, this is not ideal if sizex < CUDA_BLOCK_SIZE
    unsigned nblock = (sizex + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
    struct dim3 dg{nblock, unsigned(sizey), unsigned(sizez)};
		struct dim3 db{CUDA_BLOCK_SIZE, 1, 1};
		cuda_run_kernel_4<<<dg, db>>>(sizex, sizey, sizez, sizew, kernel);

		assert(cudaGetLastError() == CUDA_SUCCESS);
		
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
						 [itlist = begin(list)] __device__ (auto ii){
							 atomicAdd((unsigned long long int*) &(itlist[ii]), (unsigned long long int) ii + 1);
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
						 [itlist = begin(list)] __device__ (auto ii, auto jj){
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][0]), (unsigned long long int) ii + 1);
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][1]), (unsigned long long int) jj + 1);
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

	size_t check_run(size_t size1, size_t size2, size_t size3){

		math::array<size_t, 4> list({size1, size2, size3, 3}, 0);

		gpu::run(size1, size2, size3,
						 [itlist = begin(list)] __device__ (auto ii, auto jj, auto kk){
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][0]), (unsigned long long int) ii + 1);
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][1]), (unsigned long long int) jj + 1);
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][2]), (unsigned long long int) kk + 1);
						 });
		
		size_t diff = 0;
		for(size_t ii = 0; ii < size1; ii++) {
			for(size_t jj = 0; jj < size2; jj++) {
				for(size_t kk = 0; kk < size3; kk++) {
					diff += ii + 1 - list[ii][jj][kk][0];
					diff += jj + 1 - list[ii][jj][kk][1];
					diff += kk + 1 - list[ii][jj][kk][2];
				}
			}
		}

		return diff;

	}
	
	size_t check_run(size_t size1, size_t size2, size_t size3, size_t size4){

		math::array<size_t, 5> list({size1, size2, size3, size4, 4}, 0);

		gpu::run(size1, size2, size3, size4,
						 [itlist = begin(list)] __device__ (auto ii, auto jj, auto kk, auto ll){
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][ll][0]), (unsigned long long int) ii + 1);
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][ll][1]), (unsigned long long int) jj + 1);
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][ll][2]), (unsigned long long int) kk + 1);
							 atomicAdd((unsigned long long int*) &(itlist[ii][jj][kk][ll][3]), (unsigned long long int) ll + 1);
						 });
		
		size_t diff = 0;
		for(size_t ii = 0; ii < size1; ii++) {
			for(size_t jj = 0; jj < size2; jj++) {
				for(size_t kk = 0; kk < size3; kk++) {
					for(size_t ll = 0; ll < size4; ll++) {
						diff += ii + 1 - list[ii][jj][kk][ll][0];
						diff += jj + 1 - list[ii][jj][kk][ll][1];
						diff += kk + 1 - list[ii][jj][kk][ll][2];
						diff += ll + 1 - list[ii][jj][kk][ll][3];
					}
				}
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

	SECTION("3D"){
		REQUIRE(gpu::check_run(2, 2, 2) == 0);
		REQUIRE(gpu::check_run(7, 2, 2) == 0);
		REQUIRE(gpu::check_run(7, 57, 57) == 0);
		REQUIRE(gpu::check_run(32, 23, 18) == 0);
		REQUIRE(gpu::check_run(213, 27, 78) == 0);
		REQUIRE(gpu::check_run(2500, 10, 12) == 0);
		REQUIRE(gpu::check_run(7, 1023, 12) == 0);	
		REQUIRE(gpu::check_run(1, 11, 1229) == 0);	
	}
	
	SECTION("4D"){
		REQUIRE(gpu::check_run(2, 2, 2, 2) == 0);
		REQUIRE(gpu::check_run(7, 2, 2, 2) == 0);
		REQUIRE(gpu::check_run(7, 57, 57, 57) == 0);
		REQUIRE(gpu::check_run(32, 23, 45, 18) == 0);
		REQUIRE(gpu::check_run(35, 213, 27, 78) == 0);
		REQUIRE(gpu::check_run(2500, 10, 11, 12) == 0);
		REQUIRE(gpu::check_run(7, 1023, 11, 12) == 0);	
		REQUIRE(gpu::check_run(1, 1, 11, 1229) == 0);	
	}
}

#endif
#endif
#endif
