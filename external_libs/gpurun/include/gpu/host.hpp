/* -*- indent-tabs-mode: t -*- */

#ifndef GPURUN__GPU__HOST
#define GPURUN__GPU__HOST

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#ifdef ENABLE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

#include<cstddef> // std::size_t
#include<cassert>
#include <iostream>

#ifdef ENABLE_GPU
#define GPU_FUNCTION __host__ __device__
#define GPU_LAMBDA __device__
#else
#define GPU_FUNCTION
#define GPU_LAMBDA
#endif

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>

namespace gpu {

#ifdef ENABLE_GPU
template <class ErrorType>
void check_error(ErrorType const & error){
#ifdef ENABLE_CUDA
	if(error != cudaError_t(CUDA_SUCCESS)){
		std::cout << "**************************************************************************\n\n";
		std::cout << "  CUDA ERROR: '" << cudaGetErrorString(error) << "'.\n";
		std::cout << "\n**************************************************************************\n" << std::endl;		
		abort();
	}
#endif
#ifdef ENABLE_HIP
	if(error != hipError_t(hipSuccess)){
		std::cout << "**************************************************************************\n\n";
		std::cout << "  HIP ERROR: '" << hipGetErrorString(error) << "'.\n";
		std::cout << "\n**************************************************************************\n" << std::endl;
		abort();
	}
#endif
}
#endif

auto id() {
#ifdef ENABLE_GPU
#ifdef ENABLE_CUDA
	cudaDeviceProp prop;
	check_error(cudaGetDeviceProperties(&prop, 0));
	auto uuid = prop.uuid;
#endif
#ifdef ENABLE_HIP
	hipUUID uuid;
	check_error(hipDeviceGetUuid(&uuid, 0));
#endif
	using namespace boost::archive::iterators;
	using it = base64_from_binary<transform_width<unsigned char*, 6, 8>>;
	return std::string(it((unsigned char*)&uuid), it((unsigned char*)&uuid + sizeof(uuid)));
#else
	return 0;
#endif
}

auto last_error() {
#ifdef ENABLE_CUDA
	return cudaGetLastError();
#endif
#ifdef ENABLE_HIP
	return hipGetLastError();
#endif
}

template <typename KernelType>
auto max_blocksize(KernelType const & kernel){
	[[maybe_unused]] int mingridsize = 0;
	int blocksize = 0;
#ifdef ENABLE_CUDA
	check_error(cudaOccupancyMaxPotentialBlockSize(&mingridsize, &blocksize, kernel));
#endif
#ifdef ENABLE_HIP
	check_error(hipOccupancyMaxPotentialBlockSize(&mingridsize, &blocksize, kernel));
#endif
	return blocksize;
}

void sync(){
#ifdef ENABLE_CUDA
	check_error(cudaStreamSynchronize(0));
#endif
#ifdef ENABLE_HIP
	check_error(hipStreamSynchronize(0));
#endif
}

auto get_current_device() {
	int device = 0;
#ifdef ENABLE_CUDA
	check_error(cudaGetDevice(&device));
#endif
	return device;
}

template <typename PointerType>
auto get_device(PointerType pointer) {
	int device = 0;
#ifdef ENABLE_CUDA
	cudaPointerAttributes attr{};
	check_error(cudaPointerGetAttributes(&attr, raw_pointer_cast(pointer)));
	assert(attr.type == cudaMemoryTypeManaged);
	device = attr.device;
#endif
	return device;
}

constexpr auto cpu_device() {
#ifdef ENABLE_CUDA
	return cudaCpuDeviceId;
#endif
	return 0;
}

template <typename PointerType, typename SizeType>
void prefetch_to_device(PointerType pointer, SizeType byte_count, int device) {
#ifdef ENABLE_CUDA
	check_error(cudaMemPrefetchAsync(raw_pointer_cast(pointer), byte_count, device));
#endif
}

}
#endif

#ifdef GPURUN__HOST__UNIT_TEST
#undef GPURUN__HOST__UNIT_TEST

TEST_CASE(GPURUN_TEST_FILE, GPURUN_TEST_TAG) {

}

#endif
