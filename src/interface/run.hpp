/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__RUN
#define INQ__INTERFACE__RUN

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

namespace inq {
namespace interface {

struct {

	std::string name() const {
		return "run";
	}

	std::string one_line() const {
		return "Runs the simulation.";
	}
	
	auto ground_state() const{
		auto ions = systems::ions::load(".default_ions");
		systems::electrons electrons(ions, options::electrons::load(".default_electrons_options"));
		
		if(not electrons.try_load(".default_orbitals")){
			ground_state::initial_guess(ions, electrons);
		}
		auto result = ground_state::calculate(ions, electrons, options::theory::load(".default_theory"));
		electrons.save(".default_orbitals");
		return result;
	}
} const run;

}
}
#endif

#ifdef INQ_INTERFACE_RUN_UNIT_TEST
#undef INQ_INTERFACE_RUN_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
