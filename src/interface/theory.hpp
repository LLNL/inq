/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__THEORY
#define INQ__INTERFACE__THEORY

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
		return "theory";
	}

	std::string one_line() const {
		return "Defines the theory used to represent the electrons-electron interaction.";
	}
	
	void non_interacting() const{
		auto theo = options::theory{}.non_interacting();
		theo.save(input::environment::global().comm(), ".default_theory");
	}
} const theory;

}
}
#endif

#ifdef INQ_INTERFACE_THEORY_UNIT_TEST
#undef INQ_INTERFACE_THEORY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
