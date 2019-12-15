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

TEST_CASE("Test kripton", "[test::kripton]") {

	using namespace Catch::literals;

	std::vector<input::atom> geo;

	geo.push_back( "Kr" | math::d3vector(0.0, 0.0, 0.0));
		
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);

	SECTION("Non interacting"){
		
		input::config conf;

		conf.extra_states = 3;

		systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), input::interaction::non_interacting(), conf);
		
		auto energy = electrons.calculate_ground_state();
		
		/*
			OCTOPUS RESULTS: (Spacing 0.350877)

			Eigenvalues [H]
			#st  Spin   Eigenvalue      Occupation
			1   --    -4.509167       2.000000
			2   --    -3.852590       2.000000
			3   --    -3.852590       2.000000
			4   --    -3.852590       2.000000
			5   --    -2.864000       0.000000
			6   --    -2.863818       0.000000
			7   --    -2.863600       0.000000
			
			Energy [H]:
      Total       =       -32.13387478
      Free        =       -32.13387478
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =       -32.13387478
      Hartree     =         0.00000000
      Int[n*v_xc] =         0.00000000
      Exchange    =         0.00000000
      Correlation =         0.00000000
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         0.00000000
      -TS         =        -0.00000000
      Kinetic     =         9.31957760
      External    =       -41.45345238
      Non-local   =         4.84381034

		*/
		
		//THE ION-ION ENERGY IS WRONG, IT SHOULD BE ZERO

		REQUIRE(energy.eigenvalues     == -61.740948878692_a);
		REQUIRE(energy.total()         == -82.051773872674_a);
		REQUIRE(energy.external        ==  -3.603296578514_a);
		REQUIRE(energy.nonlocal        == -17.878544893202_a);
		REQUIRE(energy.kinetic()       ==  35.547251972126_a);
		
	}

}
