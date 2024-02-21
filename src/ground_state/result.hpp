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

		utils::create_directory(comm, dirname);
		
		comm.barrier();

    energy.save(comm, dirname + "/energy");
    
		auto exception_happened = true;
		if(comm.root()) {
			
			utils::save_value(comm, dirname + "/total_iter",    total_iter,    error_message);
			utils::save_value(comm, dirname + "/dipole",        dipole,        error_message);
      utils::save_value(comm, dirname + "/magnetization", magnetization, error_message);
			utils::save_value(comm, dirname + "/num_atoms",     forces.size(), error_message);
			utils::save_array(comm, dirname + "/forces",        forces,        error_message);
			
			exception_happened = false;
			comm.broadcast_value(exception_happened);
			
		} else {
			comm.broadcast_value(exception_happened);
			if(exception_happened) throw std::runtime_error(error_message);
		}
		
		comm.barrier();
	}
	
  static auto load(std::string const & dirname) {
    auto error_message = "INQ error: Cannot load the energy from directory '" + dirname + "'.";

    result res;
    res.energy = energy_type::load(dirname + "/energy");
    
    utils::load_value(dirname + "/total_iter",     res.total_iter,     error_message);
    utils::load_value(dirname + "/dipole",         res.dipole,         error_message);
    utils::load_value(dirname + "/magnetization",  res.magnetization,  error_message);

    int num_atoms;
    utils::load_value(dirname + "/num_atoms",      num_atoms,          error_message);

    res.forces = ground_state::result::forces_type(num_atoms);
    utils::load_array(dirname + "/forces",         res.forces,         error_message);
    
    return res;
	}

  template<class OStream>
  friend OStream & operator<<(OStream & out, result const & self){

    using namespace magnitude;
    
    std::cout << "Ground-state result:\n";
    std::cout << " iterations     = " << self.total_iter << '\n';
    std::cout << " dipole         = " << self.dipole << '\n';
    std::cout << " magnetization  = " << self.magnetization << '\n';
    std::cout << " total energy   = "
              << utils::num_to_str("%.8f", self.energy.total()) << " Ha | "
              << utils::num_to_str("%.8f", self.energy.total()/in_atomic_units(1.0_eV)) << " eV \n";
    std::cout << std::endl;
    return out;
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
  res.total_iter    = 333;
  res.dipole        = vector3<double>{1.55, 2.55, 3.55};
  res.magnetization = vector3<double>{-1.55, -2.55, -3.55};
  res.energy.ion(1.55);
	res.energy.ion_kinetic(2.55);
	res.energy.eigenvalues(3.55);
	res.energy.external(4.55);
	res.energy.non_local(5.55);
	res.energy.hartree(6.55);
	res.energy.xc(7.55);
	res.energy.nvxc(8.55);
	res.energy.exact_exchange(10.55);
  res.forces = ground_state::result::forces_type{65, vector3<double>{3.55, 4.55, 5.55}};

  CHECK(res.total_iter == 333);
  CHECK(res.dipole  == vector3<double>{1.55, 2.55, 3.55});
  CHECK(res.magnetization == vector3<double>{-1.55, -2.55, -3.55});
	CHECK(res.energy.ion() == 1.55);
	CHECK(res.energy.ion_kinetic() == 2.55);
	CHECK(res.energy.eigenvalues() == 3.55);
	CHECK(res.energy.external() == 4.55);
	CHECK(res.energy.non_local() == 5.55);
	CHECK(res.energy.hartree() == 6.55);
	CHECK(res.energy.xc() == 7.55);
	CHECK(res.energy.nvxc() == 8.55);
	CHECK(res.energy.exact_exchange() == 10.55);
  CHECK(res.forces == ground_state::result::forces_type{65, vector3<double>{3.55, 4.55, 5.55}});
  
  res.save(comm, "result_save");
  auto read_res = ground_state::result::load("result_save");

  CHECK(read_res.total_iter == 333);
  CHECK(read_res.dipole  == vector3<double>{1.55, 2.55, 3.55});
  CHECK(read_res.magnetization == vector3<double>{-1.55, -2.55, -3.55});
	CHECK(read_res.energy.ion() == 1.55);
	CHECK(read_res.energy.ion_kinetic() == 2.55);
	CHECK(read_res.energy.eigenvalues() == 3.55);
	CHECK(read_res.energy.external() == 4.55);
	CHECK(read_res.energy.non_local() == 5.55);
	CHECK(read_res.energy.hartree() == 6.55);
	CHECK(read_res.energy.xc() == 7.55);
	CHECK(read_res.energy.nvxc() == 8.55);
	CHECK(read_res.energy.exact_exchange() == 10.55);
  CHECK(read_res.forces == ground_state::result::forces_type{65, vector3<double>{3.55, 4.55, 5.55}});

}
#endif

