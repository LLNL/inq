/* -*- indent-tabs-mode: t -*- */

#ifndef OPTIONS__REAL_TIME
#define OPTIONS__REAL_TIME

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <magnitude/time.hpp>
#include <utils/merge_optional.hpp>

#include <optional>
#include <cassert>


namespace inq {
namespace options {

class real_time {

public:

	enum class electron_propagator { ETRS, CRANK_NICOLSON };

private:

	std::optional<double> dt_;
	std::optional<int> num_steps_;
	std::optional<electron_propagator> prop_;

public:
	
	auto dt(quantity<magnitude::time> dt) const {
		real_time solver = *this;;
		solver.dt_ = dt.in_atomic_units();
		return solver;
	}

	auto dt() const {
		return dt_.value_or(0.01);
	}

	auto num_steps(double etol) const {
		real_time solver = *this;;
		solver.num_steps_ = etol;
		return solver;
	}
				
	auto num_steps() const {
		return num_steps_.value_or(100);
	}

	auto etrs() {
		real_time solver = *this;;
		solver.prop_ = electron_propagator::ETRS;
		return solver;
	}

	auto crank_nicolson() const {
		real_time solver = *this;;
		solver.prop_ = electron_propagator::CRANK_NICOLSON;
		return solver;
	}
	
	auto propagator() const {
		return prop_.value_or(electron_propagator::ETRS);
	}

};
    
}
}
#endif

#ifdef INQ_OPTIONS_REAL_TIME_UNIT_TEST
#undef INQ_OPTIONS_REAL_TIME_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
  using namespace inq::magnitude;	
  using namespace Catch::literals;

	SECTION("Defaults"){

    options::real_time solver;

    CHECK(solver.dt() == 0.01_a);
    CHECK(solver.num_steps() == 100);
    CHECK(solver.propagator() == options::real_time::electron_propagator::ETRS);		
    
  }

  SECTION("Composition"){

    auto solver = options::real_time{}.num_steps(1000).dt(0.05_atomictime).crank_nicolson();
    
    CHECK(solver.num_steps() == 1000);
    CHECK(solver.dt() == 0.05_a);
		CHECK(solver.propagator() == options::real_time::electron_propagator::CRANK_NICOLSON);
  }

}
#endif
