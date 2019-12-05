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

namespace utils {

  template <class comm_type>
  class distribution {

  public:

		distribution(const long size, const comm_type & comm):
      size_(size),
      comm_(comm)
    {
      auto bsize = (size_ + comm_.size() - 1)/comm_.size();

      start_ = bsize*comm_.rank();
      end_ = std::min(bsize*(comm_.rank() + 1), size_);

      assert(local_size() <= bsize);
		}

    auto size() const {
      return size_;
    }

    auto local_size() const {
      return end_ - start_;
    }
    
	protected:

    long size_;
    comm_type comm_;
    long start_;
    long end_;
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("class utils::distribution", "[utils::distribution]") {
  
  using namespace Catch::literals;
  using math::d3vector;

  const int NN = 1000;

  auto comm = boost::mpi3::environment::get_world_instance();
  
  utils::distribution<boost::mpi3::communicator> dist(NN, comm);

  REQUIRE(NN == dist.size());

  auto calculated_size = dist.local_size();

  comm.all_reduce_in_place_n(&calculated_size, 1, std::plus<>{});

  REQUIRE(NN == calculated_size);
  
}
#endif

    
#endif
