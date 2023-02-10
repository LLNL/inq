/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__ARRAY_ITERATOR
#define INQ__PARALLEL__ARRAY_ITERATOR

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


#include <parallel/partition.hpp>
#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/detail/datatype.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

template <typename ArrayType>
class array_iterator {

  partition part_;
  mutable parallel::cartesian_communicator<1> comm_;
  ArrayType arr_;
  long istep_;

  struct end_type {
  };
  
public:

  array_iterator(partition part, parallel::cartesian_communicator<1> comm, ArrayType const & arr):
    part_(std::move(part)),
    comm_(std::move(comm)),
    arr_(part.block_size()),
    istep_(0)
  {
		assert(arr.size() == part_.local_size());
		arr_({0, part_.local_size()}) = arr;
  }

  auto operator!=(end_type) const {
    return istep_ != comm_.size();
  }

  void operator++(){
    auto next_proc = comm_.rank() + 1;
    if(next_proc == comm_.size()) next_proc = 0;
    auto prev_proc = comm_.rank() - 1;
    if(prev_proc == -1) prev_proc = comm_.size() - 1;

    auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();

    if(istep_ < comm_.size() - 1){
      MPI_Sendrecv_replace(raw_pointer_cast(arr_.data_elements()), arr_.num_elements(), mpi_type, prev_proc, istep_, next_proc, istep_, comm_.get(), MPI_STATUS_IGNORE);
    }
    
    istep_++;
  }

  auto ipart() const {
    auto ip = istep_ + comm_.rank();
    if(ip >= comm_.size()) ip -= comm_.size();
    return ip;
  }

  auto operator*() const {
    return arr_({0, part_.local_size(ipart())});
  }

  auto operator->() const {
    return &arr_({0, part_.local_size(ipart())});
  }
  
  static auto end() {
    return end_type{};
  }
  
};
}
}
#endif

#ifdef INQ_PARALLEL_ARRAY_ITERATOR_UNIT_TEST
#undef INQ_PARALLEL_ARRAY_ITERATOR_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>
#include <math/array.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class parallel::array_iterator", "[parallel::array_iterator]") {
  
	using namespace inq;
	using namespace Catch::literals;
  
	auto comm = parallel::cartesian_communicator<1>{boost::mpi3::environment::get_world_instance(), {}};

  long size = 12345;
  
  parallel::partition part(size, comm);

  math::array<double, 1> arr(part.local_size(), double(comm.rank() + 1.0));

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
#endif
