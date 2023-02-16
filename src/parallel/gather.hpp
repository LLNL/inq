/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__GATHER
#define INQ__PARALLEL__GATHER

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

#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/detail/datatype.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

template <class ArrayType, class PartType, class CommType>
auto gather(ArrayType const & array, PartType const & part, CommType & comm, int root) {

  if(comm.size() == 1) {
    return array;
  } else {

    ArrayType ret;
    
    auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();
    
    std::vector<int> recvcounts(comm.size());
    std::vector<int> displs(comm.size());
    
    if(comm.rank() == root){
      ret.reextent(part.size());
      for(int ipart = 0; ipart < comm.size(); ipart++){
        recvcounts[ipart] = part.local_size(ipart);
        displs[ipart] = part.start(ipart);
      }
    }
		
    MPI_Gatherv(raw_pointer_cast(array.data_elements()), part.local_size(), mpi_type, raw_pointer_cast(ret.data_elements()), recvcounts.data(), displs.data(), mpi_type, root, comm.get());
    
    return ret;
  }
}

}
}

#endif

#ifdef INQ_PARALLEL_GATHER_UNIT_TEST
#undef INQ_PARALLEL_GATHER_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>
#include <math/array.hpp>

#include <mpi3/environment.hpp>
#include <parallel/arbitrary_partition.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
  
	auto comm = parallel::cartesian_communicator<1>{boost::mpi3::environment::get_world_instance(), {}};

  auto nn = 15;
  parallel::arbitrary_partition part(nn*(comm.rank() + 1), comm);
  
  math::array<double, 1> local_array(part.local_size(), double(comm.rank() + 1.0));

  auto array = gather(local_array, part, comm, comm.size() - 1);

  if(comm.rank() == comm.size() - 1){
    int index = 0;
    for(int ipart = 0; ipart < comm.size(); ipart++){
      for(int ii = 0; ii < nn*(ipart + 1); ii++){
        CHECK(array[index] == ipart + 1.0);
        index++;
      }
    }
  }
  
}
#endif
