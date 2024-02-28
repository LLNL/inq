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

  std::string dirname_;

public:

  long total_steps;
  double total_time;
	std::vector<double> time;
	std::vector<double> total_energy;
	
  results(std::string const & arg_dirname):
    dirname_(arg_dirname),
    total_steps(0),
    total_time(0.0){
  }

  template <class ObservablesType>
  void operator()(ObservablesType const & observables){
    
    total_steps = observables.iter() + 1;
    total_time = observables.time();
    time.push_back(observables.time());
		total_energy.push_back(observables.energy().total());

  }

	void save(parallel::communicator & comm) const {
		auto error_message = "INQ error: Cannot save real_time::results to directory '" + dirname_ + "'.";

    utils::create_directory(comm, dirname_);
		utils::save_value(comm, dirname_ + "/total_steps",    total_steps,    error_message);
		utils::save_value(comm, dirname_ + "/total_time",     total_time,     error_message);
		utils::save_array(comm, dirname_ + "/time",           time,           error_message);
		utils::save_array(comm, dirname_ + "/total_energy",   total_energy,   error_message);
	}
  
  static auto load(std::string const & dirname) {
    auto error_message = "INQ error: Cannot load real_time::results from directory '" + dirname + "'.";

    results res(dirname);

    utils::load_value(dirname + "/total_steps",     res.total_steps,     error_message);
    utils::load_value(dirname + "/total_time",      res.total_time,      error_message);
		
		res.time.resize(res.total_steps + 1);
		utils::load_array(dirname + "/time",            res.time,            error_message);

		res.total_energy.resize(res.total_steps + 1);
		utils::load_array(dirname + "/total_energy",    res.total_energy,    error_message);
		
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

