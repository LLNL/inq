/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__PARTITION
#define INQ__UTILS__PARTITION

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


#include <utils/global_index.hpp>

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace utils {

class partition {
	
public:

	auto local_size() const {
		return end_ - start_;
	}
	
	partition(const long size, int comm_size = 1, int comm_rank = 0)
		:comm_size_(comm_size),
		 size_(size)
	{
		
		bsize_ = (size_ + comm_size_ - 1)/comm_size_;
		
		if(size_ > 0) assert(bsize_ > 0);

		start_ = std::min(bsize_*comm_rank, size_);
		end_ = std::min(bsize_*(comm_rank + 1), size_);
		
		assert(local_size() <= bsize_);
		assert(end_ >= start_);
		assert(end_ <= size);
	}
	
	partition(const long size, const boost::mpi3::communicator & comm)
		:partition(size, comm.size(), comm.rank()){
	}
	
	auto operator*=(const long factor) {
		size_ *= factor;
		start_ *= factor;
		end_ *= factor;
		bsize_ *= factor;
		
		return *this;
	}

	friend auto operator*(const long factor, partition part){
		part *= factor;
		return part;
	}
	
	auto size() const {
		return size_;
	}
	
	auto start() const {
		return start_;
	}

	auto start(int part) const {
		return std::min(bsize_*part, size_);
	}
	
	auto end() const {
		return end_;
	}
	
	auto end(int part) const {
		return std::min(bsize_*(part + 1), size_);
	}
	
	auto local_size(int part) const {
		return end(part) - start(part);
	}
	
	auto parallel() const {
		return comm_size_ > 1;
	}
	
	auto contains(long index) const {
		return start() <= index and index < end();
	}
	
	constexpr auto local_to_global(long local_i) const {
		return global_index(start_ + local_i);
	}
	
	constexpr auto global_to_local(global_index global_i) const {
		return global_i.value() - start_;
	}
	
	auto comm_size() const {
		return comm_size_;
	}

	auto block_size() const {
		return bsize_;
	}
	
	auto location(long global_i) const {
		return global_i/bsize_;
	}
	
protected:
	
	long comm_size_;
	long size_;
	long start_;
	long end_;
	long bsize_;
  
};
}
}

#ifdef INQ_UTILS_PARTITION_UNIT_TEST
#undef INQ_UTILS_PARTITION_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unitcell.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class utils::partition", "[utils::partition]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;

  const int NN = 1033;

  auto comm = boost::mpi3::environment::get_world_instance();
  
	inq::utils::partition part(NN, comm);

  auto next = comm.rank() + 1;
  if(next == comm.size()) next = 0;
  
  auto prev = comm.rank() - 1;
  if(prev == -1) prev = comm.size() - 1;

  SECTION("Total"){

    CHECK(NN == part.size());

    auto calculated_size = part.local_size();
    
    comm.all_reduce_in_place_n(&calculated_size, 1, std::plus<>{});
    
    CHECK(NN == calculated_size);
  }

  SECTION("Upper bound"){
    auto boundary_value = part.end();
    
    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ next, /* source = */ prev, 0, 0);
    
    if(comm.rank() != 0){
      CHECK(boundary_value == part.start());
    } else {
      CHECK(boundary_value == NN);
    }
  }

  SECTION("Lower bound"){
    auto boundary_value = part.start();

    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ prev, /* source = */ next, 1, 1);
    
    if(comm.rank() != comm.size() - 1){
      CHECK(boundary_value == part.end());
    } else {
      CHECK(boundary_value == 0);
    }
  }

	SECTION("Location"){

		for(long ig = part.start(); ig < part.end(); ig++){
			CHECK(part.location(ig) == comm.rank());
		}
		
	}

	long factor = 13;

	part *= factor;
	
	SECTION("Scaled - Total"){
		
    CHECK(NN*factor == part.size());
		
    auto calculated_size = part.local_size();
    
    comm.all_reduce_in_place_n(&calculated_size, 1, std::plus<>{});
    
    CHECK(NN*factor == calculated_size);
  }
	

  SECTION("Scaled - Upper bound"){
    auto boundary_value = part.end();
    
    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ next, /* source = */ prev, 0, 0);
    
    if(comm.rank() != 0){
      CHECK(boundary_value == part.start());
    } else {
      CHECK(boundary_value == NN*factor);
    }
  }

  SECTION("Scaled - Lower bound"){
    auto boundary_value = part.start();

    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ prev, /* source = */ next, 1, 1);
    
    if(comm.rank() != comm.size() - 1){
      CHECK(boundary_value == part.end());
    } else {
      CHECK(boundary_value == 0);
    }
  }
	
	SECTION("Scaled - Location"){
			
		for(long ig = part.start(); ig < part.end(); ig++){
			CHECK(part.location(ig) == comm.rank());
		}
		
	}

	SECTION("Small sizes"){
		for(int ii = 0; ii < 20; ii++) inq::utils::partition part(ii, comm);
	}

	SECTION("Check partition sizes 16 in 4"){
		inq::utils::partition part(16, 4, 0);

		CHECK(part.block_size() == 4);
		CHECK(part.local_size() == 4);
		
		CHECK(part.local_size(0) == 4);
		CHECK(part.local_size(1) == 4);
		CHECK(part.local_size(2) == 4);
		CHECK(part.local_size(3) == 4);

		CHECK(part.start(0) == 0);
		CHECK(part.start(1) == 4);
		CHECK(part.start(2) == 8);
		CHECK(part.start(3) == 12);

		CHECK(part.end(0) == 4);
		CHECK(part.end(1) == 8);
		CHECK(part.end(2) == 12);
		CHECK(part.end(3) == 16);		
		
	}
	
	SECTION("Check partition sizes 16 in 5"){
		inq::utils::partition part(16, 5, 0);

		CHECK(part.block_size() == 4);
		CHECK(part.local_size() == 4);
		
		CHECK(part.local_size(0) == 4);
		CHECK(part.local_size(1) == 4);
		CHECK(part.local_size(2) == 4);
		CHECK(part.local_size(3) == 4);
		CHECK(part.local_size(4) == 0);

		CHECK(part.start(0) == 0);
		CHECK(part.start(1) == 4);
		CHECK(part.start(2) == 8);
		CHECK(part.start(3) == 12);
		CHECK(part.start(4) == 16);

		CHECK(part.end(0) == 4);
		CHECK(part.end(1) == 8);
		CHECK(part.end(2) == 12);
		CHECK(part.end(3) == 16);
		CHECK(part.end(4) == 16);
		
	}

	SECTION("Check partition sizes 17 in 5"){
		inq::utils::partition part(17, 5, 0);

		CHECK(part.block_size() == 4);
		CHECK(part.local_size() == 4);
		
		CHECK(part.local_size(0) == 4);
		CHECK(part.local_size(1) == 4);
		CHECK(part.local_size(2) == 4);
		CHECK(part.local_size(3) == 4);
		CHECK(part.local_size(4) == 1);

		CHECK(part.start(0) == 0);
		CHECK(part.start(1) == 4);
		CHECK(part.start(2) == 8);
		CHECK(part.start(3) == 12);
		CHECK(part.start(4) == 16);

		CHECK(part.end(0) == 4);
		CHECK(part.end(1) == 8);
		CHECK(part.end(2) == 12);
		CHECK(part.end(3) == 16);
		CHECK(part.end(4) == 17);
		
	}

}
#endif

    
#endif
