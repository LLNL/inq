/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__ENVIRONMENT
#define INQ__INPUT__ENVIRONMENT

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <input/distribution.hpp>
#include <utils/merge_optional.hpp>
#include <utils/profiling.hpp>

#include <mpi3/environment.hpp>

#include <optional>
#include <cassert>

namespace inq {
namespace input {

class environment {

  public:

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
			
			if(not threaded() and base_comm_.rank() == 0){
				calimgr_.flush(); // write performance results
			}
		}

		auto dist() const {
			return distribution(base_comm_);
		}
		
  private:
		
    boost::mpi3::environment mpi_env_;
		cali::ConfigManager calimgr_;
    mutable boost::mpi3::communicator base_comm_;		

};

}
}

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_ENVIRONMENT_UNIT_TEST
#undef INQ_INPUT_ENVIRONMENT_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("input::environment", "[input::environment]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
