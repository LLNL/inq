/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__REAL_TIME__RESULTS
#define INQ__REAL_TIME__RESULTS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/vector3.hpp>
#include <hamiltonian/energy.hpp>

namespace inq {
namespace real_time {

class results {

public:

	options::real_time::observables_type obs;
  long total_steps;
  double total_time;
	std::vector<double> time;
	std::vector<double> total_energy;
	std::vector<vector3<double>> dipole;
	std::vector<vector3<double>> current;
	
  results():
    total_steps(0),
    total_time(0.0){
  }

  template <class ObservablesType>
  void operator()(ObservablesType const & observables){

    total_steps = observables.iter() + 1;
    total_time = observables.time();
    time.push_back(observables.time());
		total_energy.push_back(observables.energy().total());

		if(obs.find(options::real_time::observables::dipole) != obs.end()){
			dipole.emplace_back(observables.dipole());
		}

		if(obs.find(options::real_time::observables::current) != obs.end()){
			current.emplace_back(observables.current());
		}

		assert(total_steps == (long) time.size());
		
		if(not observables.every(500)) return;

		save(observables.electrons().full_comm(), ".inq/default_checkpoint/observables");
		observables.electrons().save(".inq/default_checkpoint/orbitals");
  }

	template <typename Comm>
	void save(Comm & comm, std::string const & dirname) const {
		auto error_message = "INQ error: Cannot save real_time::results to directory '" + dirname + "'.";

		utils::create_directory(comm, dirname);
		utils::save_value(comm, dirname + "/total_steps",    total_steps,    error_message);
		utils::save_value(comm, dirname + "/total_time",     total_time,     error_message);
		utils::save_container(comm, dirname + "/time",           time,           error_message);
		utils::save_container(comm, dirname + "/total_energy",   total_energy,   error_message);
		utils::save_container(comm, dirname + "/dipole",         dipole,         error_message);
		utils::save_container(comm, dirname + "/current",        current,        error_message);
	}
  
  static auto load(std::string const & dirname) {
    auto error_message = "INQ error: Cannot load real_time::results from directory '" + dirname + "'.";

    results res;

    utils::load_value(dirname + "/total_steps",     res.total_steps,     error_message);
    utils::load_value(dirname + "/total_time",      res.total_time,      error_message);
		
		res.time.resize(res.total_steps);
		utils::load_array(dirname + "/time",            res.time,            error_message);
		assert(res.time[res.time.size() - 2] < res.time[res.time.size() - 1]);
		assert((long) res.time.size() == res.total_steps or res.time.size() == 0ul);

		res.total_energy.resize(res.total_steps);
		utils::load_array(dirname + "/total_energy",    res.total_energy,    error_message);
		assert((long) res.total_energy.size() == res.total_steps or res.total_energy.size() == 0ul);

		utils::load_vector(dirname + "/dipole",         res.dipole);
		assert((long) res.dipole.size() == res.total_steps or res.dipole.size() == 0ul);

		utils::load_vector(dirname + "/current",        res.current);
		assert((long) res.current.size() == res.total_steps or res.current.size() == 0ul);

    return res;
	}

  template<class OStream>
  friend OStream & operator<<(OStream & out, results const & self){

    using namespace magnitude;
    
    std::cout << "Real-time results:\n";
    std::cout << "  total steps          = " << self.total_steps << '\n';
    std::cout << "  simulated time       = " << self.total_time << " atu | " <<  self.total_time/in_atomic_units(1.0_fs) << " fs \n";    
    std::cout << std::endl;
    return out;
  }
  
};

}
}
#endif

#ifdef INQ_REAL_TIME_RESULTS_UNIT_TEST
#undef INQ_REAL_TIME_RESULTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

}

#endif

