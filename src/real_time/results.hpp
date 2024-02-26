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
  long total_steps_;

public:

  results(std::string const & arg_dirname):
    dirname_(arg_dirname),
    total_steps_(0){
  }

  template <class ObservablesType>
  void operator()(ObservablesType const & observables){

    total_steps_ = observables.iter();
    
  }

	void save(parallel::communicator & comm) const {
		auto error_message = "INQ error: Cannot save the real_time::results to directory '" + dirname_ + "'.";
    
    utils::create_directory(comm, dirname_);
		utils::save_value(comm, dirname_ + "/total_steps",    total_steps_,    error_message);
    
	}
  
  static auto load(std::string const & dirname) {
    auto error_message = "INQ error: Cannot load the energy from directory '" + dirname + "'.";

    results res(dirname);

    utils::load_value(dirname + "/total_steps",     res.total_steps_,     error_message);
    
    return res;
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

