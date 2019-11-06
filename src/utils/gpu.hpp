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

#include <cuda.h>

#define CUDA_BLOCK_SIZE 1024

namespace gpu {

#ifdef HAVE_CUDA
  template <class kernel_type>
  __global__ void cuda_run_kernel(long size, kernel_type kernel){
    auto ii = blockIdx.x*blockDim.x + threadIdx.x;
    if(ii < size) kernel(ii);
  }
#endif
  
  template <class kernel_type>
  void run(long size, kernel_type kernel){
    
    auto nblock = size/CUDA_BLOCK_SIZE;
    if(nblock*CUDA_BLOCK_SIZE != size) nblock++;
    
#ifdef HAVE_CUDA
    cuda_run_kernel<<<nblock, CUDA_BLOCK_SIZE>>>(size, kernel);
    
    cudaDeviceSynchronize();
#endif
    
  }

}

#endif
