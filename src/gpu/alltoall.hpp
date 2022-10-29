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

#include <cstdlib>

#ifdef ENABLE_NCCL
#define ncclRemoteError 1347895789
#include <mpi3/nccl/communicator.hpp>
#endif

#include <math/array.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <inq_config.h>
#include <mpi.h>
#include <parallel/communicator.hpp>


namespace inq {
namespace gpu {

template <typename ArrayType>
void alltoall(ArrayType & buf, parallel::communicator & comm){
	CALI_CXX_MARK_FUNCTION;

	using type = typename ArrayType::element_type;
	auto mpi_type = boost::mpi3::detail::basic_datatype<type>();
	auto count = buf[0].num_elements();

	auto method = std::getenv("INQ_COMM");

	if(method == NULL or method == std::string("collective")){

		{
			CALI_CXX_MARK_SCOPE("alltoall:mpi");
			MPI_Alltoall(MPI_IN_PLACE, count, mpi_type, raw_pointer_cast(buf.data_elements()), count, mpi_type, comm.get());
		}
		
	} else if(method == std::string("point")) {

#ifdef HAVE_MPI_ISENDRECV_REPLACE

		std::vector<MPI_Request> reqs(comm.size(), MPI_REQUEST_NULL);

		CALI_CXX_MARK_SCOPE("alltoall:mpi");
		
		for(int iproc = 0; iproc < comm.size(); iproc++){
			MPI_Isendrecv_replace(raw_pointer_cast(buf[iproc].base()), count, mpi_type, iproc, iproc, iproc, comm.rank(), comm.get(), &reqs[iproc]);
		}
		
#else

		ArrayType copy(buf);
		std::vector<MPI_Request> reqs(comm.size()*2, MPI_REQUEST_NULL);
		
		CALI_CXX_MARK_SCOPE("alltoall:mpi");
		
		for(int iproc = 0; iproc < comm.size(); iproc++){
			MPI_Irecv(raw_pointer_cast(buf[iproc].base()), count, mpi_type, iproc, comm.rank(), comm.get(), &reqs[2*iproc]);
			MPI_Isend(raw_pointer_cast(copy[iproc].base()), count, mpi_type, iproc, iproc, comm.get(), &reqs[2*iproc + 1]);
		}

#endif

		std::vector<MPI_Status> stats(reqs.size());
		MPI_Waitall(reqs.size(), reqs.data(), stats.data());
		
	} else if(method == std::string("nccl")) {

#ifndef ENABLE_NCCL
		assert(false and "inq was compiled without nccl support");		
#else
		boost::mpi3::nccl::communicator ncomm{comm};
		ArrayType copy(buf);

		CALI_CXX_MARK_SCOPE("alltoall:nccl");

		ncclGroupStart();
		for(int iproc = 0; iproc < comm.size(); iproc++){
			ncclRecv(raw_pointer_cast(buf[iproc].base()), count*sizeof(type)/sizeof(double), ncclDouble, iproc, &ncomm, 0);
			ncclSend(raw_pointer_cast(copy[iproc].base()), count*sizeof(type)/sizeof(double), ncclDouble, iproc, &ncomm, 0);
		}
		ncclGroupEnd();

#endif

	} else {
		assert(false and "uknown communication method");		
	}	
}

}
}

#ifdef INQ_GPU_ALLTOALL_UNIT_TEST
#undef INQ_GPU_ALLTOALL_UNIT_TEST

#include <gpu/run.hpp>
#include <mpi3/environment.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("function gpu::alltoall", "[gpu::alltoall]"){

  using namespace inq;
  using namespace Catch::literals;
  
  auto comm = boost::mpi3::environment::get_world_instance();
  
  int blocksize = 100;
  
  math::array<int, 2> buffer({comm.size(), blocksize}, comm.rank());

  gpu::alltoall(buffer, comm);
	gpu::sync();

  for(int iproc = 0; iproc < comm.size(); iproc++){
    for(int ib = 0; ib < blocksize; ib++){
			CHECK(buffer[iproc][ib] == iproc);
    }
  }

}


#endif
#endif
