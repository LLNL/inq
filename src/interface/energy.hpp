/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__ENERGY
#define INQ__INTERFACE__ENERGY

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <hamiltonian/energy.hpp>

namespace inq {
namespace interface {

struct {		

	std::string name() const {
		return "energy";
	}

	std::string one_line() const {
		return "Get information about the energy obtained by a ground-state calculation.";
	}

	void operator()() const {
		auto ener = hamiltonian::energy::load(".inq/default_energy");
		std::cout << ener;
	}
	
	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {

		if(args.size() == 0) {
			operator()();
			exit(0);
    }			

	}
	
} const energy;

}
}
#endif

#ifdef INQ_INTERFACE_ENERGY_UNIT_TEST
#undef INQ_INTERFACE_ENERGY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
