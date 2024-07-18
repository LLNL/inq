/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__BLOCK_ARRAY_ITERATOR
#define INQ__PARALLEL__BLOCK_ARRAY_ITERATOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <parallel/partition.hpp>
#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/detail/datatype.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

template <typename ArrayType, typename PartitionType>
class block_array_iterator {

  PartitionType part_;
  mutable parallel::cartesian_communicator<1> comm_;
	gpu::array<typename ArrayType::element_type, 2> arr_;
  long istep_;
	long bsize_;
	
  struct end_type {
  };
  
public:

  block_array_iterator(long block_size, PartitionType part, parallel::cartesian_communicator<1> comm, ArrayType const & arr):
    part_(std::move(part)),
    comm_(std::move(comm)),
    arr_({block_size, part.max_local_size()}),
    istep_(0),
		bsize_(block_size)                  
  {
		assert(arr.size() == block_size);
		assert(arr.rotated().size() == part_.local_size());
		arr_({0, bsize_}, {0, part_.local_size()}) = arr;
  }

  auto operator!=(end_type) const {
    return istep_ != comm_.size();
  }

  void operator++(){
    auto next_proc = comm_.rank() + 1;
    if(next_proc == comm_.size()) next_proc = 0;
    auto prev_proc = comm_.rank() - 1;
    if(prev_proc == -1) prev_proc = comm_.size() - 1;

		if constexpr(not is_vector3<typename ArrayType::element_type>::value){
			auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();
			
			if(istep_ < comm_.size() - 1){
				MPI_Sendrecv_replace(raw_pointer_cast(arr_.data_elements()), arr_.num_elements(), mpi_type, prev_proc, istep_, next_proc, istep_, comm_.get(), MPI_STATUS_IGNORE);
			}
		} else {
			using base_type = typename ArrayType::element_type::element_type;
			auto mpi_type = boost::mpi3::detail::basic_datatype<base_type>();
			
			if(istep_ < comm_.size() - 1){
				MPI_Sendrecv_replace((base_type *) raw_pointer_cast(arr_.data_elements()), 3*arr_.num_elements(), mpi_type, prev_proc, istep_, next_proc, istep_, comm_.get(), MPI_STATUS_IGNORE);
			}
		}
			
    istep_++;
  }

  auto ipart() const {
    auto ip = istep_ + comm_.rank();
    if(ip >= comm_.size()) ip -= comm_.size();
    return ip;
  }

  auto operator*() const {
    return arr_({0, bsize_}, {0, part_.local_size(ipart())});
  }

  auto operator->() const {
    return &arr_({0, bsize_}, {0, part_.local_size(ipart())});
  }
  
  static auto end() {
    return end_type{};
  }
  
};
}
}
#endif

#ifdef INQ_PARALLEL_BLOCK_ARRAY_ITERATOR_UNIT_TEST
#undef INQ_PARALLEL_BLOCK_ARRAY_ITERATOR_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <gpu/array.hpp>

#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
  
	auto comm = parallel::cartesian_communicator<1>{boost::mpi3::environment::get_world_instance(), {}};

	long bsize = 12;
  long size = 12345;
  
  parallel::partition part(size, comm);

	SECTION("double"){
		
		gpu::array<double, 2> arr({bsize, part.local_size()}, double(comm.rank() + 1.0));
		
		auto ipart = comm.rank();
		for(parallel::block_array_iterator pai(bsize, part, comm, arr); pai != pai.end(); ++pai){

			CHECK(ipart == pai.ipart());
			CHECK((*pai).size() == bsize);
			CHECK((*pai).rotated().size() == part.local_size(ipart));
			
			for(int ii = 0; ii < (*pai).size(); ii++){
				for(int jj = 0; jj < bsize; jj++) CHECK((*pai)[jj][ii] == ipart + 1.0);
			}
			
			ipart++;
			if(ipart == comm.size()) ipart = 0;
		}

	}

	SECTION("vector3"){
		
		gpu::array<vector3<double>, 2> arr({bsize, part.local_size()}, {comm.rank() + 1.0, comm.rank()*2.0, comm.rank() - 5.0});
		
		auto ipart = comm.rank();
		for(parallel::block_array_iterator pai(bsize, part, comm, arr); pai != pai.end(); ++pai){
			
			CHECK(ipart == pai.ipart());
			CHECK( (*pai).size() == bsize);
			CHECK( (*pai).rotated().size() == part.local_size(ipart));
			
			for(int ii = 0; ii < (*pai).size(); ii++){
				for(int jj = 0; jj < bsize; jj++) {
					CHECK((*pai)[jj][ii][0] == ipart + 1.0);
					CHECK((*pai)[jj][ii][1] == ipart*2.0);
					CHECK((*pai)[jj][ii][2] == ipart - 5.0);
				}
			}
			
			ipart++;
			if(ipart == comm.size()) ipart = 0;
		}

	}
    
}
#endif
