/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__CELL
#define INQ__INTERFACE__CELL

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
		return "cell";
	}

	std::string one_line() const {
		return "Defines the simulation cell.";
	}
	
	void operator()(){
		auto cell = systems::ions::load(".default_ions").cell();
		if(input::environment::global().comm().root()) std::cout << cell;
	}
	
	static void cubic(quantity<magnitude::length> const aa, int periodicity = 3){
		systems::ions ions(systems::cell::cubic(aa).periodicity(periodicity));
		ions.save(input::environment::global().comm(), ".default_ions");
	}

} cell ;

}
}
#endif

#ifdef INQ_INTERFACE_CELL_UNIT_TEST
#undef INQ_INTERFACE_CELL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
