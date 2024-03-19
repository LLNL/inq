/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__ARRAY_ITERATOR
#define INQ__PARALLEL__ARRAY_ITERATOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <parallel/block_array_iterator.hpp>

namespace inq{
namespace parallel {

template <typename ArrayType, typename PartitionType>
class array_iterator {

	using base_iterator_type = block_array_iterator<gpu::array<typename ArrayType::element_type, 2>, PartitionType>;

	base_iterator_type bit_;
  
public:

  array_iterator(PartitionType part, parallel::cartesian_communicator<1> comm, ArrayType const & arr):
		bit_(1, part, comm, arr.partitioned(1)){		
  }

	template <typename EndType>
  auto operator!=(EndType const & end) const {
    return bit_ != end;
  }

  void operator++(){
		++bit_;
  }

  auto ipart() const {
		return bit_.ipart();
  }

  auto operator*() const {
    return bit_->flatted();
  }

  auto operator->() const {
    return &(bit_->flatted());
  }
  
  static auto end() {
    return base_iterator_type::end();
  }
  
};
}
}
#endif

#ifdef INQ_PARALLEL_ARRAY_ITERATOR_UNIT_TEST
#undef INQ_PARALLEL_ARRAY_ITERATOR_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <gpu/array.hpp>

#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
  
	auto comm = parallel::cartesian_communicator<1>{boost::mpi3::environment::get_world_instance(), {}};

  long size = 12345;
  
  parallel::partition part(size, comm);

	SECTION("double"){
		
		gpu::array<double, 1> arr(part.local_size(), double(comm.rank() + 1.0));
		
		auto ipart = comm.rank();
		for(parallel::array_iterator pai(part, comm, arr); pai != pai.end(); ++pai){
			
			CHECK(ipart == pai.ipart());
			CHECK(pai->size() == part.local_size(ipart)); 
			
			for(int ii = 0; ii < pai->size(); ii++){
				CHECK((*pai)[ii] == ipart + 1.0);
			}
			
			ipart++;
			if(ipart == comm.size()) ipart = 0;
		}

	}

	SECTION("vector3"){
		
		gpu::array<vector3<double>, 1> arr(part.local_size(), {comm.rank() + 1.0, comm.rank()*2.0, comm.rank() - 5.0});
		
		auto ipart = comm.rank();
		for(parallel::array_iterator pai(part, comm, arr); pai != pai.end(); ++pai){
			
			CHECK(ipart == pai.ipart());
			CHECK(pai->size() == part.local_size(ipart)); 
			
			for(int ii = 0; ii < pai->size(); ii++){
				CHECK((*pai)[ii][0] == ipart + 1.0);
				CHECK((*pai)[ii][1] == ipart*2.0);
				CHECK((*pai)[ii][2] == ipart - 5.0);
			}
			
			ipart++;
			if(ipart == comm.size()) ipart = 0;
		}

	}
    
}
#endif
