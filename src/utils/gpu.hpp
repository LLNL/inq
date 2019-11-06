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

#ifdef HAVE_CUDA
#include <cuda.h>
#endif

#define CUDA_BLOCK_SIZE 1024

namespace gpu {

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
    unsigned nblock = (size + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
    
    cuda_run_kernel_1<<<nblock, CUDA_BLOCK_SIZE>>>(size, kernel);
    
    cudaDeviceSynchronize();
#endif
    
  }

#ifdef HAVE_CUDA
  template <class kernel_type>
  __global__ void cuda_run_kernel_2(unsigned sizex, unsigned sizey, kernel_type kernel){
    auto ix = blockIdx.x*blockDim.x + threadIdx.x;
    auto iy = blockIdx.y*blockDim.y + threadIdx.y;
    if(ix < sizex && iy < sizey) kernel(ix, iy);
  }
#endif
 
  template <class kernel_type>
  void run(long sizex, long sizey, kernel_type kernel){

#ifdef HAVE_CUDA
    //OPTIMIZATION, this is not ideal if sizex < CUDA_BLOCK_SIZE
    unsigned nblock = (sizex + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
    
    cuda_run_kernel_2<<<{nblock, unsigned(sizey)}, {CUDA_BLOCK_SIZE, 1}>>>(sizex, sizey, kernel);
    
    cudaDeviceSynchronize();
#endif
    
  }

}

#endif
