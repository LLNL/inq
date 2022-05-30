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

  boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {1, boost::mpi3::fill});

	auto basis_comm = cart_comm.axis(1);

  for(auto ecut = 25.0_Ha; ecut <= 300.0_Ha; ecut += 25.0_Ha){
    
    systems::box box = systems::box::cubic(6.66_b).cutoff_energy(ecut);
    basis::real_space rs(box, basis_comm);
    
    for(int ist = 1; ist <= 256; ist *= 2){
      
      basis::field_set<basis::real_space, complex> phi(rs, ist, cart_comm);
      
      auto iter_start_time = std::chrono::high_resolution_clock::now();
      
      for(int irep = 0; irep < reps; irep++){
        auto fphi = operations::space::to_fourier(phi);
      }
      
      auto new_time = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_seconds = new_time - iter_start_time;
    
      std::cout << ecut.in_atomic_units() << '\t' << ist << '\t' << elapsed_seconds.count()/reps*1000.0 << "\tms" << std::endl;

    }
  }
  
}
