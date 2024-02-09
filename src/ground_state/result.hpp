/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__RESULT
#define INQ__GROUND_STATE__RESULT

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/vector3.hpp>
#include <hamiltonian/energy.hpp>

namespace inq {
namespace ground_state {

struct result {

  using energy_type = hamiltonian::energy;
  using forces_type = gpu::array<vector3<double>, 1>;

  int total_iter;
  vector3<double> dipole;
  vector3<double> magnetization;
  energy_type energy;
  forces_type forces;

	void save(parallel::communicator & comm, std::string const & dirname) const {
		auto error_message = "INQ error: Cannot save the ground_state::result to directory '" + dirname + "'.";
		
		comm.barrier();

    energy.save(comm, dirname + "/energy");
    
		auto exception_happened = true;
		if(comm.root()) {
			
			try { std::filesystem::create_directories(dirname); }
			catch(...) {
				comm.broadcast_value(exception_happened);
				throw std::runtime_error(error_message);
			}
      
			utils::save_value(comm, dirname + "/dipole",        dipole,        error_message);
      utils::save_value(comm, dirname + "/magnetization", magnetization, error_message);
			utils::save_value(comm, dirname + "/total_iter",    total_iter,    error_message);
			utils::save_array(comm, dirname + "/forces",        forces,        error_message);
			
			exception_happened = false;
			comm.broadcast_value(exception_happened);
			
		} else {
			comm.broadcast_value(exception_happened);
			if(exception_happened) throw std::runtime_error(error_message);
		}
		
		comm.barrier();
	}
	
  
  
};

}
}
#endif

#ifdef INQ_GROUND_STATE_RESULT_UNIT_TEST
#undef INQ_GROUND_STATE_RESULT_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
  
  ground_state::result res;
  res.total_iter    = 666;
  res.dipole        = vector3<double>{1.0, 2.0, 3.0};
  res.magnetization = vector3<double>{-1.0, -2.0, -3.0};
  res.energy.ion(1.0);
	res.energy.ion_kinetic(2.0);
	res.energy.eigenvalues(3.0);
	res.energy.external(4.0);
	res.energy.non_local(5.0);
	res.energy.hartree(6.0);
	res.energy.xc(7.0);
	res.energy.nvxc(8.0);
	res.energy.exact_exchange(10.0);
  res.forces = ground_state::result::force_type{65, vector3<double>{3.0, 4.0, 5.0}};

  res.save(comm, "result_save");
  
}
#endif

