/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <utils/match.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

  int const reps = 10;

	auto & env = inq::input::environment::global();
	auto comm = env.comm();
	comm.nccl_init();
	
  if(comm.root()) printf("#  size [MB]      p2p [GB/s] agregate [GB/s]\n");
  
  for(int blocksize = 62; blocksize < 80*1024*1024 ; blocksize *= sqrt(2)){
    
    gpu::array<complex, 2> buffer({comm.size(), blocksize}, double(comm.rank()));

    parallel::alltoall(buffer, comm);

    auto ok = true;
    for(int iproc = 0; iproc < comm.size(); iproc++){
      for(int ib = 0; ib < blocksize; ib++){
        ok = ok and buffer[iproc][ib] == double(iproc);
      }
    }

    if(not ok){
      std::cerr << "alltoall failed for size " << blocksize << std::endl;
      exit(1);
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    for(int irep = 0; irep < reps; irep++) parallel::alltoall(buffer, comm);
    std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;

    auto ttime = comm.all_reduce_value(time.count(), boost::mpi3::plus<>{});
    ttime /= comm.size();
    
    auto rate = reps*blocksize*2*sizeof(complex)/ttime/1e9;
    
    if(comm.root()) printf("%12.4f\t%12.2f\t%12.2f\n", sizeof(complex)*blocksize/1e6, rate, rate*comm.size()*(comm.size() - 1.0)/2);
    
  }
  
}
