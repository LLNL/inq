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
	
	real_time(){
	}

	static auto dt(quantity<magnitude::time> dt) {
		real_time solver;
		solver.dt_ = dt.in_atomic_units();
		return solver;
	}

	auto dt() const {
		return dt_.value_or(0.01);
	}

	auto static num_steps(double etol) {
		real_time solver;
		solver.num_steps_ = etol;
		return solver;
	}
				
	auto num_steps() const {
		return num_steps_.value_or(100);
	}

	auto static etrs() {
		real_time solver;
		solver.prop_ = electron_propagator::ETRS;
		return solver;
	}

	auto static crank_nicolson() {
		real_time solver;
		solver.prop_ = electron_propagator::CRANK_NICOLSON;
		return solver;
	}
	
	auto propagator() const {
		return prop_.value_or(electron_propagator::ETRS);
	}
	
	friend auto operator|(const real_time & solver1, const real_time & solver2){
		using utils::merge_optional;

		real_time rsolver;
		rsolver.dt_	= merge_optional(solver1.dt_, solver2.dt_);
		rsolver.num_steps_	= merge_optional(solver1.num_steps_, solver2.num_steps_);
		rsolver.prop_	= merge_optional(solver1.prop_, solver2.prop_);		
		return rsolver;
	}
    
private:

	std::optional<double> dt_;
	std::optional<int> num_steps_;
	std::optional<electron_propagator> prop_;
		
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

    auto solver = options::real_time::num_steps(1000) | options::real_time::dt(0.05_atomictime) | options::real_time::crank_nicolson();
    
    CHECK(solver.num_steps() == 1000);
    CHECK(solver.dt() == 0.05_a);
		CHECK(solver.propagator() == options::real_time::electron_propagator::CRANK_NICOLSON);
  }

}
#endif
