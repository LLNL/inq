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

#include <cstdlib>

#include <math/array.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <inq_config.h>
#include <mpi.h>
#include <mpi3/communicator.hpp>


namespace inq {
namespace gpu {

template <typename ArrayType>
void alltoall(ArrayType & buf, boost::mpi3::communicator & comm){

	auto mpi_type = boost::mpi3::detail::basic_datatype<typename ArrayType::element_type>();
	auto count = buf[0].num_elements();

	auto method = std::getenv("INQ_COMM");

	if(method == NULL or method == std::string("collective")){
		
		MPI_Alltoall(MPI_IN_PLACE, count, mpi_type, raw_pointer_cast(buf.data_elements()), count, mpi_type, comm.get());

	} else if(method == std::string("point")) {
		
		ArrayType copy(buf);
		
		std::vector<MPI_Request> reqs(comm.size()*2, MPI_REQUEST_NULL);

		//We should use MPI_Isendrecv_replace here but it is not implemented in some OpenMPI versions 
		for(int iproc = 0; iproc < comm.size(); iproc++){
			MPI_Irecv(raw_pointer_cast(buf[iproc].base()), count, mpi_type, iproc, comm.rank(), comm.get(), &reqs[2*iproc]);
			MPI_Isend(raw_pointer_cast(copy[iproc].base()), count, mpi_type, iproc, iproc, comm.get(), &reqs[2*iproc + 1]);
		}
		
		std::vector<MPI_Status> stats(comm.size()*2);
		MPI_Waitall(reqs.size(), reqs.data(),  stats.data());

	} else {
		assert(false and "uknown communication method");		
	}	
}

}
}

#ifdef INQ_GPU_ALLTOALL_UNIT_TEST
#undef INQ_GPU_ALLTOALL_UNIT_TEST

#include <mpi3/environment.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("function gpu::alltoall", "[gpu::alltoall]"){

  using namespace inq;
  using namespace Catch::literals;
  
  auto comm = boost::mpi3::environment::get_world_instance();
  
  int blocksize = 100;
  
  math::array<int, 2> buffer({comm.size(), blocksize}, comm.rank());
  
  gpu::alltoall(buffer, comm);
  
  for(int iproc = 0; iproc < comm.size(); iproc++){
    for(int ib = 0; ib < blocksize; ib++){
      CHECK(buffer[iproc][ib] == iproc);
    }
  }

}


#endif
#endif
