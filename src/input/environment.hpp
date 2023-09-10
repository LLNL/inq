/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__ENVIRONMENT
#define INQ__INPUT__ENVIRONMENT

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/parallelization.hpp>
#include <utils/merge_optional.hpp>
#include <utils/profiling.hpp>

#include <spdlog/spdlog.h>

#include <inq_config.h>  // for ENABLE_CUDA macro

#include <mpi3/environment.hpp>

#ifdef ENABLE_CUDA
#include <multi/array.hpp>  // for multi::array
#include <thrust/system/cuda/memory.h>  // for ::thrust::cuda::allocator
#include <csetjmp>
#include <csignal>
#endif

#include <cassert>
#include <optional>

namespace inq {
namespace input {

class environment {

 private:
		static auto & threaded_impl() {
			static bool threaded_ = false;

			return threaded_;
		}

	public:

		static auto const & threaded() {
			return threaded_impl();
		}

    environment(int argc, char** argv, bool use_threads = false):
      mpi_env_(argc, argv, use_threads?boost::mpi3::thread_level::multiple:boost::mpi3::thread_level::single),
			base_comm_(mpi_env_.get_world_instance())
    {
			if(use_threads){
				assert(mpi_env_.thread_support() == boost::mpi3::thread_level::multiple);
			}
			
			threaded_impl() = use_threads;
			
			if(not use_threads and base_comm_.rank() == 0){
				calimgr_.add("runtime-report");
				calimgr_.start();
			}

			CALI_MARK_BEGIN("inq_environment");
    }

		~environment(){

			CALI_MARK_END("inq_environment");

			base_comm_.barrier();

			if(not threaded() and base_comm_.rank() == 0){
				calimgr_.flush(); // write performance results
			}

		}

		auto par() const {
			return parallelization(base_comm_);
		}

	auto world() {return mpi_env_.world();}

  private:
		
    boost::mpi3::environment mpi_env_;
		cali::ConfigManager calimgr_;
    mutable parallel::communicator base_comm_;		
};

}
}
#endif

#ifdef INQ_INPUT_ENVIRONMENT_UNIT_TEST
#undef INQ_INPUT_ENVIRONMENT_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
