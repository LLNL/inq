/* -*- indent-tabs-mode: t -*- */

#ifndef UTILS__DISTRIBUTION
#define UTILS__DISTRIBUTION

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

#include <cassert>
#include <array>

#include <mpi3/communicator.hpp>

namespace utils {

  class distribution {

  public:

    auto local_size() const {
      return end_ - start_;
    }
		
		distribution(const long size, const boost::mpi3::communicator & comm):
			comm_size_(comm.size()),
      size_(size){
			
      auto bsize = (size_ + comm_size_ - 1)/comm_size_;

			if(size_ > 0) assert(bsize > 0);

      start_ = bsize*comm.rank();
      end_ = std::min(bsize*(comm.rank() + 1), size_);

      assert(local_size() <= bsize);
			assert(end_ >= start_);
		}

		auto operator*=(const long factor) {
			size_ *= factor;
			start_ *= factor;
			end_ *= factor;
			
			return *this;
		}
		
    auto size() const {
      return size_;
    }

    auto start() const {
      return start_;
    }

    auto end() const {
      return end_;
    }
    
		auto parallel() const {
			return comm_size_ > 1;
		}

		auto contains(long index) const {
			return start() <= index and index < end();
		}

		auto local_to_global(long local_i) const {
			return start_ + local_i;
		}

		auto global_to_local(long global_i) const {
			return global_i - start_;
		}
		
	protected:

		long comm_size_;
    long size_;
    long start_;
    long end_;
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class utils::distribution", "[utils::distribution]") {
  
  using namespace Catch::literals;
  using math::vec3d;

  const int NN = 1033;

  auto comm = boost::mpi3::environment::get_world_instance();
  
  utils::distribution dist(NN, comm);

  auto next = comm.rank() + 1;
  if(next == comm.size()) next = 0;
  
  auto prev = comm.rank() - 1;
  if(prev == -1) prev = comm.size() - 1;

  SECTION("Total"){

    REQUIRE(NN == dist.size());

    auto calculated_size = dist.local_size();
    
    comm.all_reduce_in_place_n(&calculated_size, 1, std::plus<>{});
    
    REQUIRE(NN == calculated_size);
  }

  SECTION("Upper bound"){
    auto boundary_value = dist.end();
    
    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ next, /* source = */ prev, 0, 0);
    
    if(comm.rank() != 0){
      REQUIRE(boundary_value == dist.start());
    } else {
      REQUIRE(boundary_value == NN);
    }
  }

  SECTION("Lower bound"){
    auto boundary_value = dist.start();

    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ prev, /* source = */ next, 1, 1);
    
    if(comm.rank() != comm.size() - 1){
      REQUIRE(boundary_value == dist.end());
    } else {
      REQUIRE(boundary_value == 0);
    }
  }


	long factor = 13;

	dist *= factor;
	
	SECTION("Scaled - Total"){
		
    REQUIRE(NN*factor == dist.size());
		
    auto calculated_size = dist.local_size();
    
    comm.all_reduce_in_place_n(&calculated_size, 1, std::plus<>{});
    
    REQUIRE(NN*factor == calculated_size);
  }

  SECTION("Scaled - Upper bound"){
    auto boundary_value = dist.end();
    
    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ next, /* source = */ prev, 0, 0);
    
    if(comm.rank() != 0){
      REQUIRE(boundary_value == dist.start());
    } else {
      REQUIRE(boundary_value == NN*factor);
    }
  }

  SECTION("Scaled - Lower bound"){
    auto boundary_value = dist.start();

    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ prev, /* source = */ next, 1, 1);
    
    if(comm.rank() != comm.size() - 1){
      REQUIRE(boundary_value == dist.end());
    } else {
      REQUIRE(boundary_value == 0);
    }
  }
}
#endif

    
#endif
