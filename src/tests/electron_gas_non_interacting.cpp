/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <catch2/catch.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>

TEST_CASE("Test non interacting electron gas", "[test::non_interacting_electron_gas]") {

	using namespace Catch::literals;

	systems::ions ions(input::cell::cubic(10.0, 10.0, 10.0));

	input::config conf;
	conf.extra_states = 2;
	conf.excess_charge = 14.0;
	conf.theory = input::electronic_theory::NON_INTERACTING;
		
	systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), conf);

	auto energy = electrons.calculate_ground_state();

	//Octopus results are:
	// Energy: 2.36870506
	// Eigenvalues: 0.000000 0.197392 0.394784
	REQUIRE(energy.total() == 2.3687083213_a);
	REQUIRE(energy.kinetic() == 2.3687083213_a);
	REQUIRE(energy.eigenvalues == 2.3687083213_a);
	
}
