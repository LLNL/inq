/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__PARTITION
#define INQ__PARALLEL__PARTITION

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


#include <parallel/global_index.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/detail/datatype.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

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
	
	partition(const long size, const parallel::communicator & comm)
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
	
	constexpr auto start() const {
		return start_;
	}

	constexpr auto start(int part) const {
		return std::min(bsize_*part, size_);
	}
	
	constexpr auto end() const {
		return end_;
	}
	
	constexpr auto end(int part) const {
		return std::min(bsize_*(part + 1), size_);
	}
	
	constexpr auto local_size(int part) const {
		return end(part) - start(part);
	}
	
	auto parallel() const {
		return comm_size_ > 1;
	}
	
	auto contains(long index) const {
		return start() <= index and index < end();
	}

	auto contains(long index, int part) const {
		return start(part) <= index and index < end(part);
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

	auto waste() const {
		auto total_elements = block_size()*comm_size();
		return (total_elements - size())/double(size());
	}
	
	template <class ArrayType, class CommType>
	auto gather(ArrayType const & array, CommType & comm, int root) const {
		if(comm.size() == 1) {
			return array;
		} else {

			ArrayType ret;

			auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();

			std::vector<int> recvcounts(comm.size());
			std::vector<int> displs(comm.size());

			if(comm.rank() == root){
				ret.reextent(size_);
				for(int ipart = 0; ipart < comm.size(); ipart++){
					recvcounts[ipart] = local_size(ipart);
					displs[ipart] = start(ipart);
				}
			}
			
			MPI_Gatherv(raw_pointer_cast(array.data_elements()), local_size(), mpi_type, raw_pointer_cast(ret.data_elements()), recvcounts.data(), displs.data(), mpi_type, root, comm.get());

			return ret;
		}
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

#ifdef INQ_PARALLEL_PARTITION_UNIT_TEST
#undef INQ_PARALLEL_PARTITION_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class parallel::partition", "[parallel::partition]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;

  const int NN = 1033;

  auto comm = boost::mpi3::environment::get_world_instance();
  
	inq::parallel::partition part(NN, comm);

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

	SECTION("Waste"){
		if(comm.size() == 1) CHECK(part.waste() == 0.0_a);
		if(comm.size() == 2) CHECK(part.waste() == 0.0009680542_a);
		if(comm.size() == 3) CHECK(part.waste() == 0.0019361084_a);
		if(comm.size() == 4) CHECK(part.waste() == 0.0029041626_a);
		if(comm.size() == 5) CHECK(part.waste() == 0.0019361084_a);
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
		for(int ii = 0; ii < 20; ii++) inq::parallel::partition part(ii, comm);
	}

	SECTION("Check partition sizes 16 in 4"){
		inq::parallel::partition part(16, 4, 0);

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
		
		CHECK(part.waste() == 0.0_a);
	}
	
	SECTION("Check partition sizes 16 in 5"){
		inq::parallel::partition part(16, 5, 0);

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

		CHECK(part.waste() == 0.25_a);		
	}

	SECTION("Check partition sizes 17 in 5"){
		inq::parallel::partition part(17, 5, 0);

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

		CHECK(part.waste() == 0.1764705882_a);
	}

}
#endif

    
#endif
