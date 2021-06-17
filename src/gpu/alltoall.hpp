/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GPU__ALLTOALL
#define INQ__GPU__ALLTOALL

/*
 Copyright (C) 2021 Xavier Andrade

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

#include <inq_config.h>
#include <mpi.h>
#include <mpi3/communicator.hpp>

namespace inq {
namespace gpu {

void alltoall(void *recvbuf, long recvcount, MPI_Datatype recvtype, boost::mpi3::communicator & comm){
  MPI_Alltoall(MPI_IN_PLACE, recvcount, recvtype, recvbuf, recvcount, recvtype, comm.get());


  /*
			std::vector<MPI_Request> reqs(comm.size()*2);
			for(int iproc = 0; iproc < comm.size(); iproc++){
				if(iproc == comm.rank()) continue;
				MPI_Irecv(raw_pointer_cast(buffer2[iproc].base()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, iproc, comm.rank(), comm.get(), &reqs[2*iproc]);
				MPI_Isend(raw_pointer_cast(buffer[iproc].base()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, iproc, iproc, comm.get(), &reqs[2*iproc + 1]);
			}

			buffer2[comm.rank()] == buffer[comm.rank()];

			std::vector<MPI_Status> stats(comm.size()*2);
			MPI_Waitall(reqs.size(), reqs.data(),  stats.data());
  */

}

}
}

#ifdef INQ_GPU_ALLTOALL_UNIT_TEST
#undef INQ_GPU_ALLTOALL_UNIT_TEST

#include <mpi3/environment.hpp>

#include <math/array.hpp>
#include <catch2/catch.hpp>

TEST_CASE("function gpu::alltoall", "[gpu::alltoall]"){

  using namespace inq;
  using namespace Catch::literals;
  
  auto comm = boost::mpi3::environment::get_world_instance();
  
  int blocksize = 100;
  
  math::array<int, 2> buffer({comm.size(), blocksize}, comm.rank());
  
  gpu::alltoall(raw_pointer_cast(buffer.data_elements()), blocksize, MPI_INT, comm);
  
  for(int iproc = 0; iproc < comm.size(); iproc++){
    for(int ib = 0; ib < blocksize; ib++){
      CHECK(buffer[iproc][ib] == iproc);
    }
  }

}


#endif
#endif
