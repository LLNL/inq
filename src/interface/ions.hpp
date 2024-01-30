/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__IONS
#define INQ__INTERFACE__IONS

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

namespace inq {
namespace interface {

struct {
	
	std::string name() const {
		return "ions";
	}

	std::string one_line() const {
		return "Defines the ions in the simulation.";
	}

	void operator()() const {
		auto ions = systems::ions::load(".default_ions");		
		if(input::environment::global().comm().root()) std::cout << ions;
	}

	void add(input::species const & sp, vector3<quantity<magnitude::length>> const & pos) const {
		auto ions = systems::ions::load(".default_ions");
		ions.insert(sp, pos);
		ions.save(input::environment::global().comm(), ".default_ions");
	}

	void clear() const {
		auto ions = systems::ions::load(".default_ions");
		ions.clear();
		ions.save(input::environment::global().comm(), ".default_ions");
	}
	
} const ions;

}
}
#endif

#ifdef INQ_INTERFACE_IONS_UNIT_TEST
#undef INQ_INTERFACE_IONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
