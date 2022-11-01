/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2022 Xavier Andrade

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

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

  int const reps = 10;

  input::environment env(argc, argv);
	parallel::communicator comm = boost::mpi3::environment::get_world_instance();

  if(comm.root()) printf("#  size [MB]      p2p [GB/s] agregate [GB/s]\n");
  
  for(int blocksize = 62; blocksize < 80*1024*1024 ; blocksize *= sqrt(2)){
    
    math::array<complex, 2> buffer({comm.size(), blocksize}, double(comm.rank()));

    gpu::alltoall(buffer, comm);

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
    for(int irep = 0; irep < reps; irep++) gpu::alltoall(buffer, comm);
    std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;

    auto ttime = comm.all_reduce_value(time.count(), boost::mpi3::plus<>{});
    ttime /= comm.size();
    
    auto rate = reps*blocksize*2*sizeof(complex)/ttime/1e9;
    
    if(comm.root()) printf("%12.4f\t%12.2f\t%12.2f\n", sizeof(complex)*blocksize/1e6, rate, rate*comm.size()*(comm.size() - 1.0)/2);
    
  }
  
}
