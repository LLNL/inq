/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__ELECTRONS
#define INQ__INTERFACE__ELECTRONS

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
		return "electrons";
	}

	std::string one_line() const {
		return "Defines the electrons in the simulation and how they are represented.";
	}

	static void extra_states(int nstates){
		auto el_opts = options::electrons::load(".default_electrons_options").extra_states(nstates);
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}

	static void cutoff(quantity<magnitude::energy> ecut){
		auto el_opts = options::electrons::load(".default_electrons_options").cutoff(ecut);
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}

	static void fourier_pseudo(){
		auto el_opts = options::electrons::load(".default_electrons_options").fourier_pseudo();
		el_opts.save(input::environment::global().comm(), ".default_electrons_options");
	}

} electrons;

}
}
#endif

#ifdef INQ_INTERFACE_ELECTRONS_UNIT_TEST
#undef INQ_INTERFACE_ELECTRONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
