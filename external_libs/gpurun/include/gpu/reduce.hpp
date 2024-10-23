/* -*- indent-tabs-mode: t -*- */

#ifndef GPURUN__GPU__REDUCE
#define GPURUN__GPU__REDUCE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <cassert>

#ifndef ENABLE_GPU
#include <numeric>
#else
#include <thrust/execution_policy.h>
#include <thrust/transform_reduce.h>
#endif

#include <gpu/run.hpp>
#include <gpu/array.hpp>
#include <gpu/host.hpp>

namespace gpu {

struct reduce {
	explicit reduce(long arg_size):
		size(arg_size){		
	}
	long size;
};

template <typename ArrayType>
struct array_access {
  ArrayType array;

  GPU_FUNCTION auto operator()(long ii) const {
    return array[ii];
  }

  GPU_FUNCTION auto operator()(long ix, long iy) const {
    return array[ix][iy];
  }
  
};

template <typename Type, typename KernelType>
Type run(reduce const & red, Type const init, KernelType kernel) {

	auto const size = red.size;
	auto range = boost::multi::extension_t{0l, size};

#ifndef ENABLE_GPU
	return std::transform_reduce(range.begin(), range.end(), init, std::plus<>{}, kernel);
#else
	return thrust::transform_reduce(thrust::device, range.begin(), range.end(), kernel, init, std::plus<>{});
#endif
}

#ifdef ENABLE_GPU
template <typename KernelType, typename ArrayType>
__global__ void reduce_kernel_rr(long sizex, long sizey, KernelType kernel, ArrayType odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename ArrayType::element *) shared_mem;
	
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;	

	if(ix < sizex and iy < sizey){
		reduction_buffer[tid] = kernel(ix, iy);
	} else {
		reduction_buffer[tid] = (typename ArrayType::element) 0.0;
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
	if (tid == 0) odata[blockIdx.x][blockIdx.y] = reduction_buffer[0];

}
#endif

template <typename Type, typename KernelType>
Type run(gpu::reduce const & redx, gpu::reduce const & redy, Type const init, KernelType kernel) {

	auto const sizex = redx.size;	
	auto const sizey = redy.size;	
  
#ifndef ENABLE_GPU

  auto accumulator = init;
	for(long iy = 0; iy < sizey; iy++){
		for(long ix = 0; ix < sizex; ix++){
			accumulator += kernel(ix, iy);
		}
	}
  return accumulator;

#else

	const int bsizex = 1024;
	const int bsizey = 1;

	unsigned nblockx = (sizex + bsizex - 1)/bsizex;
	unsigned nblocky = (sizey + bsizey - 1)/bsizey;
	
	gpu::array<Type, 2> result({nblockx, nblocky});

	struct dim3 dg{nblockx, nblocky};
	struct dim3 db{bsizex, bsizey};

	reduce_kernel_rr<<<dg, db, bsizex*bsizey*sizeof(Type)>>>(sizex, sizey, kernel, begin(result));
  check_error(last_error());
	
  if(nblockx*nblocky == 1) {
    gpu::sync();
    return init + result[0][0];
  } else {
    return run(gpu::reduce(nblockx*nblocky), init, array_access<decltype(begin(result.flatted()))>{begin(result.flatted())});
  }
  
#endif
}

#ifdef ENABLE_GPU
template <typename KernelType, typename ArrayType>
__global__ void reduce_kernel_rrr(long sizex, long sizey, long sizez, KernelType kernel, ArrayType odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename ArrayType::element *) shared_mem;
	
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int iz = blockIdx.z*blockDim.z + threadIdx.z;	

	if(ix < sizex and iy < sizey and iz < sizez){
		reduction_buffer[tid] = kernel(ix, iy, iz);
	} else {
		reduction_buffer[tid] = (typename ArrayType::element) 0.0;
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
	if (tid == 0) odata[blockIdx.x][blockIdx.y][blockIdx.z] = reduction_buffer[0];

}
#endif

template <typename Type, typename KernelType>
Type run(reduce const & redx, reduce const & redy, reduce const & redz, Type const init, KernelType kernel) {

	auto const sizex = redx.size;	
	auto const sizey = redy.size;
	auto const sizez = redz.size;
	
	if(sizex == 0 or sizey == 0 or sizez == 0) return init;
	
#ifndef ENABLE_GPU

  auto accumulator = init;
	for(long iy = 0; iy < sizey; iy++){
		for(long ix = 0; ix < sizex; ix++){
			for(long iz = 0; iz < sizez; iz++){
				accumulator += kernel(ix, iy, iz);
			}
		}
	}
  return accumulator;
	
#else

	auto blocksize = max_blocksize(reduce_kernel_rrr<KernelType, decltype(begin(std::declval<gpu::array<Type, 3>&>()))>);
	
	const unsigned bsizex = blocksize;
	const unsigned bsizey = 1;
	const unsigned bsizez = 1;

	unsigned nblockx = (sizex + bsizex - 1)/bsizex;
	unsigned nblocky = (sizey + bsizey - 1)/bsizey;
	unsigned nblockz = (sizez + bsizez - 1)/bsizez;

	gpu::array<Type, 3> result({nblockx, nblocky, nblockz});

	struct dim3 dg{nblockx, nblocky, nblockz};
	struct dim3 db{bsizex, bsizey, bsizez};

	reduce_kernel_rrr<<<dg, db, bsizex*bsizey*bsizez*sizeof(Type)>>>(sizex, sizey, sizez, kernel, begin(result));
	check_error(last_error());

  if(nblockx*nblocky*nblockz == 1) {
    gpu::sync();
    return init + result[0][0][0];
  } else {
    return run(gpu::reduce(nblockx*nblocky*nblockz), init, array_access<decltype(begin(result.flatted().flatted()))>{begin(result.flatted().flatted())});
  }
  
#endif
}

#ifdef ENABLE_GPU
template <typename KernelType, typename ArrayType>
__global__ void reduce_kernel_vr(long sizex, long sizey, KernelType kernel, ArrayType odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename ArrayType::element *) shared_mem; // {blockDim.x, blockDim.y}
	
	// each thread loads one element from global to shared mem
  unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.y;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;

	if(ix >= sizex) return;
	
	if(iy < sizey){
		reduction_buffer[threadIdx.x + blockDim.x*tid] = kernel(ix, iy);
	} else {
		reduction_buffer[threadIdx.x + blockDim.x*tid] = (typename ArrayType::element) 0.0;
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


template <typename Type, typename KernelType>
gpu::array<Type, 1> run(long sizex, reduce const & redy, Type const init, KernelType kernel) { 

	auto const sizey = redy.size;	

#ifndef ENABLE_GPU

  gpu::array<Type, 1> accumulator(sizex, init);

  for(long iy = 0; iy < sizey; iy++){
    for(long ix = 0; ix < sizex; ix++){
      accumulator[ix] += kernel(ix, iy);
    }
  }
  
  return accumulator;
  
#else

	gpu::array<Type, 2> result;
	
	auto blocksize = max_blocksize(reduce_kernel_vr<KernelType, decltype(begin(result))>);
	
	unsigned bsizex = 4; //this seems to be the optimal value
	if(sizex <= 2) bsizex = sizex;
	unsigned bsizey = blocksize/bsizex;

	assert(bsizey > 1);
	
	unsigned nblockx = (sizex + bsizex - 1)/bsizex;
	unsigned nblocky = (sizey + bsizey - 1)/bsizey;
		
	result.reextent({nblocky, sizex});

	struct dim3 dg{nblockx, nblocky};
  struct dim3 db{bsizex, bsizey};

  auto shared_mem_size = blocksize*sizeof(Type);

  assert(shared_mem_size <= 48*1024);
  
  reduce_kernel_vr<<<dg, db, shared_mem_size>>>(sizex, sizey, kernel, begin(result));	
  check_error(last_error());
	
  if(nblocky == 1) {
    gpu::sync();

		assert(result[0].size() == sizex);

		gpu::run(result[0].size(), [res = begin(result[0]), init] GPU_LAMBDA (auto ii) {
			res[ii] += init;
		});
			
    return result[0];
  } else {
    return run(sizex, reduce(nblocky), init, array_access<decltype(begin(result.transposed()))>{begin(result.transposed())});
  }
  
#endif

}

#ifdef ENABLE_GPU
template <typename KernelType, typename ArrayType>
__global__ void reduce_kernel_vrr(long sizex, long sizey,long sizez, KernelType kernel, ArrayType odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename ArrayType::element *) shared_mem; // {blockDim.x, blockDim.y}
	
	// each thread loads one element from global to shared mem
  unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.y;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int iz = blockIdx.z*blockDim.z + threadIdx.z;

	if(ix >= sizex) return;
	
	if(iy < sizey and iz < sizez){
		reduction_buffer[threadIdx.x + blockDim.x*tid] = kernel(ix, iy, iz);
	} else {
		reduction_buffer[threadIdx.x + blockDim.x*tid] = (typename ArrayType::element) 0.0;
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
	if (tid == 0) odata[blockIdx.y][blockIdx.z][ix] = reduction_buffer[threadIdx.x];

}
#endif

template <typename Type, typename KernelType>
gpu::array<Type, 1>  run(long sizex, reduce const & redy, reduce const & redz, Type const init, KernelType kernel) {

	auto const sizey = redy.size;
	auto const sizez = redz.size;	
	
#ifndef ENABLE_GPU

  gpu::array<Type, 1> accumulator(sizex, init);

	for(long iz = 0; iz < sizez; iz++){
		for(long iy = 0; iy < sizey; iy++){
			for(long ix = 0; ix < sizex; ix++){
				accumulator[ix] += kernel(ix, iy, iz);
			}
		}
	}
  
  return accumulator;
  
#else

	gpu::array<Type, 3> result;
	
	auto blocksize = max_blocksize(reduce_kernel_vrr<KernelType, decltype(begin(result))>);
	
	unsigned bsizex = 4; //this seems to be the optimal value
	if(sizex <= 2) bsizex = sizex;
	unsigned bsizey = blocksize/bsizex;
	unsigned bsizez = 1;

	assert(bsizey > 1);
	assert(bsizex*bsizey*bsizez == blocksize);
	
	unsigned nblockx = (sizex + bsizex - 1)/bsizex;
	unsigned nblocky = (sizey + bsizey - 1)/bsizey;
	unsigned nblockz = (sizez + bsizez - 1)/bsizez;
		
	result.reextent({nblocky, nblockz, sizex});

	struct dim3 dg{nblockx, nblocky, nblockz};
  struct dim3 db{bsizex, bsizey, bsizez};

  auto shared_mem_size = blocksize*sizeof(Type);

  assert(shared_mem_size <= 48*1024);
  
  reduce_kernel_vrr<<<dg, db, shared_mem_size>>>(sizex, sizey, sizez, kernel, begin(result));	
  check_error(last_error());
	
  if(nblocky*nblockz == 1) {
    gpu::sync();

		assert(result[0][0].size() == sizex);

		gpu::run(result[0][0].size(), [res = begin(result[0][0]), init] GPU_LAMBDA (auto ii) {
			res[ii] += init;
		});
				
    return result[0][0];
  } else {
		auto && reduce_buffer = result.flatted().transposed();
    return run(sizex, reduce(nblocky*nblockz), init, array_access<decltype(begin(reduce_buffer))>{begin(reduce_buffer)});
  }
  
#endif

}

}
#endif

#ifdef GPURUN__REDUCE__UNIT_TEST
#undef GPURUN__REDUCE__UNIT_TEST

#include <mpi3/environment.hpp>
#include <catch2/catch_all.hpp>

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
  
struct prod3 {
  GPU_FUNCTION auto operator()(long ix, long iy, long iz) const {
    return double(ix)*double(iy)*double(iz);
  }
};

TEST_CASE(GPURUN_TEST_FILE, GPURUN_TEST_TAG) {
  
	using namespace Catch::literals;
	using Catch::Approx;

	SECTION("r"){
		const long maxsize = 129140163;
		
		int rank = 0;
		for(long nn = 1; nn <= maxsize; nn *= 3){
			CHECK(gpu::run(gpu::reduce(nn), -232.8, ident{}) == Approx(-232.8 + (nn*(nn - 1.0)/2.0)));
			rank++;
		}
	}

	SECTION("rr"){

		const long maxsize = 2*625;

		int rank = 0;
		for(long nx = 1; nx <= maxsize; nx *= 5){
			for(long ny = 1; ny <= maxsize; ny *= 5){

				auto res = gpu::run(gpu::reduce(nx), gpu::reduce(ny), 2.23, prod{});
				
				CHECK(typeid(decltype(res)) == typeid(double));
				CHECK(res == Approx(2.23 + nx*(nx - 1.0)/2.0*ny*(ny - 1.0)/2.0));
				rank++;
			}
		}
		
  }

	SECTION("rrr"){

		const long maxsize = 125;

		int rank = 0;
		for(long nx = 1; nx <= 10000; nx *= 10){
			for(long ny = 1; ny <= maxsize; ny *= 5){
				for(long nz = 1; nz <= maxsize; nz *= 5){
					
					auto res = gpu::run(gpu::reduce(nx), gpu::reduce(ny), gpu::reduce(nz), 17.89, prod3{});
					
					CHECK(typeid(decltype(res)) == typeid(double));
					CHECK(res == Approx(17.89 + nx*(nx - 1.0)/2.0*ny*(ny - 1.0)/2.0*nz*(nz - 1.0)/2.0));
					rank++;
				}
			}
		}
		
  }
	
	SECTION("vr"){

		const long maxsize = 390625;

		int rank = 0;
		for(long nx = 1; nx <= 10000; nx *= 10){
			for(long ny = 1; ny <= maxsize; ny *= 5){

				auto res = gpu::run(nx, gpu::reduce(ny), -7.7, [] GPU_LAMBDA (auto ix, auto iy) {return double(ix)*double(iy);});
					
				CHECK(typeid(decltype(res)) == typeid(gpu::array<double, 1>));
				CHECK(res.size() == nx);
				for(long ix = 0; ix < nx; ix++) CHECK(res[ix] == -7.7 + double(ix)*ny*(ny - 1.0)/2.0);
				rank++;
			}
		}
		
  }

	SECTION("vrr"){

		const long maxsize = 625;

		int rank = 0;
		for(long nx = 1; nx <= 10000; nx *= 10){
			for(long ny = 1; ny <= maxsize; ny *= 5){
				for(long nz = 1; nz <= maxsize; nz *= 5){
					
					auto res = gpu::run(nx, gpu::reduce(ny), gpu::reduce(nz), 10.0, prod3{});
					
					CHECK(typeid(decltype(res)) == typeid(gpu::array<double, 1>));
					
					CHECK(res.size() == nx);
					for(long ix = 0; ix < nx; ix++) CHECK(res[ix] == 10.0 + double(ix)*ny*(ny - 1.0)/2.0*nz*(nz - 1.0)/2.0);
					rank++;
				}
			}
		}
		
  }

}
#endif

