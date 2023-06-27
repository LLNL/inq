/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__ARRAY_ITERATOR_2D
#define INQ__PARALLEL__ARRAY_ITERATOR_2D

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

template <typename ArrayType, typename PartXType, typename PartYType>
class array_iterator_2d {

  PartXType partx_;
  PartYType party_;  
  mutable parallel::cartesian_communicator<2> comm_;
  ArrayType arr_;
  int step_;
  int xpart_;
  int ypart_;  
  int prev_proc_;
  int next_proc_;
  
  struct end_type {
  };
  
public:

  array_iterator_2d(PartXType partx, PartYType party, parallel::cartesian_communicator<2> comm, ArrayType const & arr):
    partx_(std::move(partx)),
    party_(std::move(party)),
    comm_(std::move(comm)),
    arr_({partx.max_local_size(), party.max_local_size()}),
    step_(0),
    xpart_(comm_.coordinates()[0]),
    ypart_(comm_.coordinates()[1])
  {

    assert(partx_.comm_size() == comm_.shape()[0]);
    assert(party_.comm_size() == comm_.shape()[1]);    
    
    arr_({0, partx_.local_size()}, {0, party_.local_size()}) = arr;

    int coords[2];

    MPI_Cart_coords(comm_.get(), comm_.rank(), 2, coords);
    coords[1]--;
    if(coords[1] == -1){
      coords[1] = comm_.shape()[1] - 1;
      coords[0]--;
      if(coords[0] == -1) coords[0] = comm_.shape()[0] - 1;
    }
    MPI_Cart_rank(comm_.get(), coords, &prev_proc_);
    
    MPI_Cart_coords(comm_.get(), comm_.rank(), 2, coords);
    coords[1]++;
    if(coords[1] == comm_.shape()[1]){
      coords[1] = 0;
      coords[0]++;
      if(coords[0] == comm_.shape()[0]) coords[0] = 0;
    }
    MPI_Cart_rank(comm_.get(), coords, &next_proc_);

  }

  auto operator!=(end_type) const {
    return step_ != comm_.size();
  }

  void operator++(){

    auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();

    if(step_ < comm_.size() - 1) {
      MPI_Sendrecv_replace(raw_pointer_cast(arr_.data_elements()), arr_.num_elements(), mpi_type, prev_proc_, step_, next_proc_, step_, comm_.get(), MPI_STATUS_IGNORE);
    }
          
    step_++;

    ypart_++;
    if(ypart_ == party_.comm_size()){
      ypart_ = 0;
      xpart_++;
      if(xpart_ == partx_.comm_size()) xpart_ = 0;
    }
  }

  auto xpart() const {
    return xpart_;
  }

  auto ypart() const {
    return ypart_;
  }

  static auto end() {
    return end_type{};
  }

  auto operator*() const {
    return arr_({0, partx_.local_size(xpart())}, {0, party_.local_size(ypart())});
  }

  auto operator->() const {
    return &arr_({0, partx_.local_size(xpart())}, {0, party_.local_size(ypart())});
  }
  
};
}
}
#endif

#ifdef INQ_PARALLEL_ARRAY_ITERATOR_2D_UNIT_TEST
#undef INQ_PARALLEL_ARRAY_ITERATOR_2D_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <gpu/array.hpp>

#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
  
	auto comm = parallel::cartesian_communicator<2>{boost::mpi3::environment::get_world_instance(), {}};

  long sizex = 1234;
  long sizey = 987;
  
  parallel::partition partx(sizex, comm.axis(0));
  parallel::partition party(sizey, comm.axis(1));

  gpu::array<double, 2> arr({partx.local_size(), party.local_size()}, comm.axis(1).rank() + 1.0 + 10000.0*(comm.axis(0).rank() + 1.0));
  
  {
    parallel::array_iterator_2d pai(partx, party, comm, arr);
    int itcount = 0;
    
    auto xpart = comm.axis(0).rank();
    for(int xstep = 0; xstep < partx.comm_size(); xstep++){
      
      auto ypart = comm.axis(1).rank();
      for(int ystep = 0; ystep < party.comm_size(); ystep++){
        
        CHECK(xpart == pai.xpart());
        CHECK(ypart == pai.ypart());

        CHECK(pai.ypart() + 1.0 + 10000.0*(pai.xpart() + 1.0) == (*pai)[0][0]);
        
        ++pai;
        ypart++;
        if(ypart == party.comm_size()){
          ypart = 0;
          xpart++;
          if(xpart == partx.comm_size()) xpart = 0;
        }

        itcount++;
      }
      
    }

    CHECK(itcount == comm.size());
  }
  
  {
    int itcount = 0;
    for(parallel::array_iterator_2d pai(partx, party, comm, arr); pai != pai.end(); ++pai){

      CHECK(pai->size() ==  partx.local_size(pai.xpart()));
      CHECK(pai->rotated().size() ==  party.local_size(pai.ypart()));
      
      for(long ix = 0; ix < partx.local_size(pai.xpart()); ix++){
        for(long iy = 0; iy < party.local_size(pai.ypart()); iy++){
          CHECK((*pai)[ix][iy] == pai.ypart() + 1.0 + 10000.0*(pai.xpart() + 1.0));
        }
      }
      
      itcount++;
    }
    
    CHECK(itcount == comm.size());
  }
  
}
#endif
