/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__ARRAY_ITERATOR_2D
#define INQ__PARALLEL__ARRAY_ITERATOR_2D

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
#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/detail/datatype.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

template <typename ArrayType>
class array_iterator_2d {

  partition partx_;
  partition party_;  
  mutable boost::mpi3::cartesian_communicator<2> comm_;
  ArrayType arr_;
  long xstep_;
  long ystep_;
  int prev_proc_;
  int next_proc_;
  
  struct end_type {
  };
  
public:

  array_iterator_2d(partition partx, partition party, boost::mpi3::cartesian_communicator<2> comm, ArrayType const & arr):
    partx_(std::move(partx)),
    party_(std::move(party)),
    comm_(std::move(comm)),
    arr_({partx.block_size(), party.block_size()}),
    xstep_(0),
    ystep_(0)      
  {
    arr_({0, partx_.local_size()}, {0, party_.local_size()}) = arr;

    int coords[2];
    MPI_Cart_coords(comm_.get(), comm_.rank(), 2, coords);
    
    assert(coords[0] == comm_.axis(0).rank());
    assert(coords[1] == comm_.axis(1).rank());

    coords[1]++;
    MPI_Cart_rank(comm_.get(), coords, &next_proc_);

    MPI_Cart_coords(comm_.get(), comm_.rank(), 2, coords);
    coords[1]--;
    MPI_Cart_rank(comm_.get(), coords, &prev_proc_);
    
  }

  auto operator!=(end_type) const {
    return xstep_ != partx_.comm_size();
  }

  void operator++(){

    auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();
    auto tag = ystep_ + xstep_*party_.comm_size();
    MPI_Sendrecv_replace(raw_pointer_cast(arr_.data_elements()), arr_.num_elements(), mpi_type, prev_proc_, tag, next_proc_, tag, comm_.get(), MPI_STATUS_IGNORE);
          
    ystep_++;
    if(ystep_ == party_.comm_size()) {
      xstep_++;
      ystep_ = 0;
    }
  }

  auto xpart() const {
    auto xp = xstep_ + comm_.coordinates()[0];
    if(xp >= partx_.comm_size()) xp -= partx_.comm_size();
    return xp;
  }

  auto ypart() const {
    auto yp = ystep_ + comm_.coordinates()[1];
    if(yp >= party_.comm_size()) yp -= party_.comm_size();
    return yp;
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

#ifdef INQ_PARALLEL_ARRAY_ITERATOR_2D_UNIT_TEST
#undef INQ_PARALLEL_ARRAY_ITERATOR_2D_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unitcell.hpp>
#include <math/array.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class parallel::array_iterator_2d", "[parallel::array_iterator_2d]") {
  
	using namespace inq;
	using namespace Catch::literals;
  
	auto comm = boost::mpi3::cartesian_communicator<2>{boost::mpi3::environment::get_world_instance(), {}};

  long sizex = 1234;
  long sizey = 987;
  
  parallel::partition partx(sizex, comm.axis(0));
  parallel::partition party(sizey, comm.axis(1));

  math::array<double, 2> arr({partx.local_size(), party.local_size()}, comm.axis(0).rank() + 1.0 + 10000.0*(comm.axis(1).rank() + 1.0));

  {
    parallel::array_iterator_2d pai(partx, party, comm, arr);
    int itcount = 0;
    
    auto xpart = comm.axis(0).rank();
    for(int xstep = 0; xstep < partx.comm_size(); xstep++){
      
      auto ypart = comm.axis(1).rank();
      for(int ystep = 0; ystep < party.comm_size(); ystep++){
        
        CHECK(xpart == pai.xpart());
        CHECK(ypart == pai.ypart());
        
        ++pai;
        
        ypart++;
        if(ypart == party.comm_size()) ypart = 0;

        itcount++;
      }
      
      xpart++;
      if(xpart == partx.comm_size()) xpart = 0;
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
          CHECK((*pai)[ix][iy] == pai.xpart() + 1.0 + 10000.0*(pai.ypart() + 1.0));
        }
      }
      
      itcount++;
    }
    
    CHECK(itcount == comm.size());
  }
  
}
#endif

    
#endif
