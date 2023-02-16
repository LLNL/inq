/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__ARBITRARY_PARTITION
#define INQ__PARALLEL__ARBITRARY_PARTITION

/*
 Copyright (C) 2023 Xavier Andrade

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

#include <parallel/global_index.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/detail/datatype.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

class arbitrary_partition {

	long local_size_;
	std::vector<long> lsizes_;
	std::vector<long> starts;
	long size_;
	long max_local_size_;
	long comm_size_;
	long start_;
	long end_;

public:
	
	auto local_size() const {
		return local_size_;
	}

	template <typename CommType>
	arbitrary_partition(long const local_size, CommType comm):
		local_size_(local_size),
		lsizes_(comm.size()),
		comm_size_(comm.size())
	{

		MPI_Allgather(&local_size, 1, MPI_LONG, lsizes_.data(), 1, MPI_LONG, comm.get());

		assert(lsizes_[comm.rank()] == local_size);
		
		size_ = 0;
		max_local_size_ = 0;
		start_ = 0;
		end_ = 0;
		
		for(unsigned ipart = 0; ipart < lsizes_.size(); ipart++){
			size_ += lsizes_[ipart];
			max_local_size_ = std::max(max_local_size_, lsizes_[ipart]);
			if(ipart < (unsigned) comm.rank()) start_ += lsizes_[ipart];
		}
		end_ = start_ + local_size_;
	}
	
	auto size() const {
		return size_;
	}
	
	constexpr auto start() const {
		return start_;
	}

	auto start(int part) const {
		long val = 0;
		for(int ipart = 0; ipart < part; ipart++) val += lsizes_[ipart];
		return val;
	}
	
	constexpr auto end() const {
		return end_;
	}

	auto end(int part) const {
		long val = 0;
		for(int ipart = 0; ipart <= part; ipart++) val += lsizes_[ipart];
		return val;
	}
	
	auto local_size(int part) const {
		return lsizes_[part];
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

	auto max_local_size() const {
		return max_local_size_;
	}
	
	auto waste() const {
		auto total_elements = max_local_size()*comm_size();
		return (total_elements - size())/double(size());
	}
	
};
}
}
#endif

#ifdef INQ_PARALLEL_ARBITRARY_PARTITION_UNIT_TEST
#undef INQ_PARALLEL_ARBITRARY_PARTITION_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
  
	parallel::communicator comm = boost::mpi3::environment::get_world_instance();
  
	auto local_size = 100 + 3*comm.rank();

	inq::parallel::arbitrary_partition part(local_size, comm);

  auto next = comm.rank() + 1;
  if(next == comm.size()) next = 0;
  
  auto prev = comm.rank() - 1;
  if(prev == -1) prev = comm.size() - 1;

  SECTION("Total"){
		
		auto total_size = local_size;

		comm.all_reduce_in_place_n(&total_size, 1, std::plus<>{});
    CHECK(total_size == part.size());

    auto calculated_size = part.local_size();
    
    comm.all_reduce_in_place_n(&calculated_size, 1, std::plus<>{});
    
    CHECK(total_size == calculated_size);

		CHECK(part.start() == part.start(comm.rank()));
		CHECK(part.end() == part.end(comm.rank()));
		
  }

	SECTION("Waste"){
		if(comm.size() == 1) CHECK(part.waste() == 0.0_a);
		if(comm.size() == 2) CHECK(part.waste() == 0.0147783251_a);
		if(comm.size() == 3) CHECK(part.waste() == 0.0291262136_a);
		if(comm.size() == 4) CHECK(part.waste() == 0.043062201_a);
		if(comm.size() == 5) CHECK(part.waste() == 0.0566037736_a);
	}

  SECTION("Upper bound"){
    auto boundary_value = part.end();
    
    comm.send_receive_replace_n(&boundary_value, 1, /* dest = */ next, /* source = */ prev, 0, 0);
    
    if(comm.rank() != 0){
      CHECK(boundary_value == part.start());
    } else {
      CHECK(boundary_value == part.size());
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

}
#endif
