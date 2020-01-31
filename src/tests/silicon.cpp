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
#include <config/path.hpp>
#include <input/atom.hpp>

TEST_CASE("Test hydrogen local pseudopotential", "[test::hydrogen_local]") {

	using namespace Catch::literals;

	std::vector<input::atom> geo;

	double a = 10.18;
	
	geo.push_back( "Si" | a*math::vec3d(0.0,  0.0,  0.0 ));
	geo.push_back( "Si" | a*math::vec3d(0.25, 0.25, 0.25));
	geo.push_back( "Si" | a*math::vec3d(0.5,  0.5,  0.0 ));
	geo.push_back( "Si" | a*math::vec3d(0.75, 0.75, 0.25));
	geo.push_back( "Si" | a*math::vec3d(0.5,  0.0,  0.5 ));
	geo.push_back( "Si" | a*math::vec3d(0.75, 0.25, 0.75));
	geo.push_back( "Si" | a*math::vec3d(0.0,  0.5,  0.5 ));
	geo.push_back( "Si" | a*math::vec3d(0.25, 0.75, 0.75));

	systems::ions ions(input::cell::cubic(a) | input::cell::finite(), geo);

	SECTION("LDA"){
		
		input::config conf;

		conf.extra_states = 4;

		systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), conf);
		
		auto energy = electrons.calculate_ground_state(input::interaction::dft());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)
		*/

		REQUIRE(energy.ion             == -0.232294220410_a);
		REQUIRE(energy.eigenvalues     == -0.234329528903_a);
		REQUIRE(energy.xc              == -0.232294220410_a);
		REQUIRE(energy.nvxc            == -0.302713349819_a);
		REQUIRE(energy.total()         == -0.446253846698_a);
		REQUIRE(energy.external        == -0.108660738870_a);
		REQUIRE(energy.nonlocal        == -0.058633055438_a);
		REQUIRE(energy.kinetic()       ==  0.416973236003_a);
		
	}

}
