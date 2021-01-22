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

struct reduce {
	explicit reduce(long arg_size):
		size(arg_size){		
	}
	long size;
};


#ifdef ENABLE_CUDA
template <class kernel_type, class array_type>
__global__ void reduce_kernel_r(long size, kernel_type kernel, array_type odata) {

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
auto run(reduce const & red, kernel_type kernel) -> decltype(kernel(0)) {

	auto const size = red.size;
	
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

  reduce_kernel_r<<<nblock, blocksize, blocksize*sizeof(type)>>>(size, kernel, begin(result));	
  check_error(cudaGetLastError());
	
  if(nblock == 1) {
    cudaDeviceSynchronize();
    return result[0];
  } else {
    return run(gpu::reduce(nblock), array_access<decltype(begin(result))>{begin(result)});
  }
  
#endif
}

#ifdef ENABLE_CUDA
template <class kernel_type, class array_type>
__global__ void reduce_kernel_rr(long sizex, long sizey, kernel_type kernel, array_type odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename array_type::element *) shared_mem;
	
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;	

	if(ix < sizex and iy < sizey){
		reduction_buffer[tid] = kernel(ix, iy);
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
	if (tid == 0) odata[blockIdx.x][blockIdx.y] = reduction_buffer[0];

}
#endif

template <class kernel_type>
auto run(reduce const & redx, reduce const & redy, kernel_type kernel) -> decltype(kernel(0, 0)) {

	auto const sizex = redx.size;	
	auto const sizey = redy.size;	
	
  using type = decltype(kernel(0, 0));
  
#ifndef ENABLE_CUDA

  type accumulator = 0.0;
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
	
	math::array<type, 2> result({nblockx, nblocky});

  reduce_kernel_rr<<<{nblockx, nblocky}, {bsizex, bsizey}, bsizex*bsizey*sizeof(type)>>>(sizex, sizey, kernel, begin(result));	
  check_error(cudaGetLastError());
	
  if(nblockx*nblocky == 1) {
    cudaDeviceSynchronize();
    return result[0][0];
  } else {
    return run(gpu::reduce(nblockx*nblocky), array_access<decltype(begin(result.flatted()))>{begin(result.flatted())});
  }
  
#endif
}

#ifdef ENABLE_CUDA
template <class kernel_type, class array_type>
__global__ void reduce_kernel_vr(long sizex, long sizey, kernel_type kernel, array_type odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename array_type::element *) shared_mem; // {blockDim.x, blockDim.y}
	
	// each thread loads one element from global to shared mem
  unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.y;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;

	if(ix >= sizex) return;
	
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
auto run(long sizex, reduce const & redy, kernel_type kernel) -> math::array<decltype(kernel(0, 0)), 1> {

	auto const sizey = redy.size;	
	
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

	math::array<type, 2> result;
	
	int mingridsize = 0;
	int blocksize = 0;

	check_error(cudaOccupancyMaxPotentialBlockSize(&mingridsize, &blocksize, reduce_kernel_vr<kernel_type, decltype(begin(result))>));
	
	unsigned bsizex = 4; //this seems to be the optimal value
	if(sizex <= 2) bsizex = sizex;
	unsigned bsizey = blocksize/bsizex;

	assert(bsizey > 1);
	
	unsigned nblockx = (sizex + bsizex - 1)/bsizex;
	unsigned nblocky = (sizey + bsizey - 1)/bsizey;
		
	result.reextent({nblocky, sizex});

	struct dim3 dg{nblockx, nblocky};
  struct dim3 db{bsizex, bsizey};

  auto shared_mem_size = blocksize*sizeof(type);

  assert(shared_mem_size <= 48*1024);
  
  reduce_kernel_vr<<<dg, db, shared_mem_size>>>(sizex, sizey, kernel, begin(result));	
  check_error(cudaGetLastError());
	
  if(nblocky == 1) {
    cudaDeviceSynchronize();

		assert(result[0].size() == sizex);
		
    return result[0];
  } else {
    return run(sizex, reduce(nblocky), array_access<decltype(begin(result.transposed()))>{begin(result.transposed())});
  }
  
#endif

}

#ifdef ENABLE_CUDA
template <class kernel_type, class array_type>
__global__ void reduce_kernel_vrr(long sizex, long sizey,long sizez, kernel_type kernel, array_type odata) {

	extern __shared__ char shared_mem[];
	auto reduction_buffer = (typename array_type::element *) shared_mem; // {blockDim.x, blockDim.y}
	
	// each thread loads one element from global to shared mem
  unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.y;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int iz = blockIdx.z*blockDim.z + threadIdx.z;

	if(ix >= sizex) return;
	
	if(iy < sizey and iz < sizez){
		reduction_buffer[threadIdx.x + blockDim.x*tid] = kernel(ix, iy, iz);
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
	if (tid == 0) odata[blockIdx.y][blockIdx.z][ix] = reduction_buffer[threadIdx.x];

}
#endif

template <class kernel_type>
auto run(long sizex, reduce const & redy, reduce const & redz, kernel_type kernel) -> math::array<decltype(kernel(0, 0, 0)), 1> {

	auto const sizey = redy.size;
	auto const sizez = redz.size;	
	
  using type = decltype(kernel(0, 0, 0));

#ifndef ENABLE_CUDA

  math::array<type, 1> accumulator(sizex, 0.0);

	for(long iz = 0; iz < sizez; iz++){
		for(long iy = 0; iy < sizey; iy++){
			for(long ix = 0; ix < sizex; ix++){
				accumulator[ix] += kernel(ix, iy, iz);
			}
		}
	}
  
  return accumulator;
  
#else

	math::array<type, 3> result;
	
	int mingridsize = 0;
	int blocksize = 0;

	check_error(cudaOccupancyMaxPotentialBlockSize(&mingridsize, &blocksize, reduce_kernel_vrr<kernel_type, decltype(begin(result))>));
	
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

  auto shared_mem_size = blocksize*sizeof(type);

  assert(shared_mem_size <= 48*1024);
  
  reduce_kernel_vrr<<<dg, db, shared_mem_size>>>(sizex, sizey, sizez, kernel, begin(result));	
  check_error(cudaGetLastError());
	
  if(nblocky*nblockz == 1) {
    cudaDeviceSynchronize();

		assert(result[0][0].size() == sizex);
		
    return result[0][0];
  } else {
		auto && reduce_buffer = result.flatted().transposed();
    return run(sizex, reduce(nblocky*nblockz), array_access<decltype(begin(reduce_buffer))>{begin(reduce_buffer)});
  }
  
#endif

}

}
}

#ifdef INQ_GPU_REDUCE_UNIT_TEST
#undef INQ_GPU_REDUCE_UNIT_TEST

#include <mpi3/environment.hpp>

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
	
struct prod3 {
  GPU_FUNCTION auto operator()(long ix, long iy, long iz) const {
    return double(ix)*double(iy)*double(iz);
  }
};

TEST_CASE("function gpu::reduce", "[gpu::reduce]") {
  
	using namespace inq;
	using namespace Catch::literals;

	auto comm = boost::mpi3::environment::get_world_instance();

	SECTION("r"){
		const long maxsize = 129140163;

		int rank = 0;
		for(long nn = 1; nn <= maxsize; nn *= 3){
			if(comm.rank() == rank%comm.size()) CHECK(gpu::run(gpu::reduce(nn), ident{}) == (nn*(nn - 1.0)/2.0));
			rank++;
		}
	}

	SECTION("rr"){

		const long maxsize = 2*625;

		int rank = 0;
		for(long nx = 1; nx <= maxsize; nx *= 5){
			for(long ny = 1; ny <= maxsize; ny *= 5){

				if(comm.rank() == rank%comm.size()){
					auto res = gpu::run(gpu::reduce(nx), gpu::reduce(ny), prod{});
					
					CHECK(typeid(decltype(res)) == typeid(double));
					CHECK(res == nx*(nx - 1.0)/2.0*ny*(ny - 1.0)/2.0);
				}
				
				rank++;
			}
		}
		
  }
	
	SECTION("vr"){

		const long maxsize = 390625;

		int rank = 0;
		for(long nx = 1; nx <= 10000; nx *= 10){
			for(long ny = 1; ny <= maxsize; ny *= 5){

				if(comm.rank() == rank%comm.size()){
					auto res = gpu::run(nx, gpu::reduce(ny), prod{});
					
					CHECK(typeid(decltype(res)) == typeid(math::array<double, 1>));
					
					CHECK(res.size() == nx);
					for(long ix = 0; ix < nx; ix++) CHECK(res[ix] == double(ix)*ny*(ny - 1.0)/2.0);
				}
				
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
					
					if(comm.rank() == rank%comm.size()){
						auto res = gpu::run(nx, gpu::reduce(ny), gpu::reduce(nz), prod3{});
						
						CHECK(typeid(decltype(res)) == typeid(math::array<double, 1>));
						
						CHECK(res.size() == nx);
						for(long ix = 0; ix < nx; ix++) CHECK(res[ix] == double(ix)*ny*(ny - 1.0)/2.0*nz*(nz - 1.0)/2.0);
					}
					
					rank++;
				}
			}
		}
		
  }

	
}

#endif
#endif

