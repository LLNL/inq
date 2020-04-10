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
	
	geo.push_back("H" | math::vec3d(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);
#if 0
	SECTION("Non interacting"){
		
		input::config conf;
		
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), conf);
		
		auto energy = electrons.calculate_ground_state(input::interaction::non_interacting());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)

			Eigenvalues [H]
			#st  Spin   Eigenvalue      Occupation
			1   --    -0.500411       1.000000
			
			Energy [H]:
      Total       =        -0.50041053
      Free        =        -0.50041053
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =        -0.50041053
      Hartree     =         0.00000000
      Int[n*v_xc] =         0.00000000
      Exchange    =         0.00000000
      Correlation =         0.00000000
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.49215618
      External    =        -0.99256671
      Non-local   =        -0.06871050

		*/
		REQUIRE(energy.eigenvalues     == -0.500213154797_a);
		REQUIRE(fabs(energy.xc)        <=  1e-10);
		REQUIRE(fabs(energy.nvxc)      <=  1e-10);
		REQUIRE(energy.total()         == -0.500213154797_a);
		REQUIRE(energy.external        == -0.125329604081_a);
		REQUIRE(energy.nonlocal        == -0.068186046673_a);
		REQUIRE(energy.kinetic()       ==  0.487958292140_a);

	}
#endif
#if 1
	SECTION("LDA"){
		
		input::config conf;
	
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), conf);
		
		auto energy = electrons.calculate_ground_state(input::interaction::dft(), input::scf_solver::density_mixing());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)
			Eigenvalues [H]
			#st  Spin   Eigenvalue      Occupation
			1   --    -0.233880       1.000000
			
			Energy [H]:
      Total      =        -0.44618708
      Free        =        -0.44618708
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =        -0.23387966
      Hartree     =         0.28283087
      Int[n*v_xc] =        -0.30315816
      Exchange    =        -0.19299390
      Correlation =        -0.03964083
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.41796657
      External    =        -0.91434707
      Non-local   =        -0.05876129


		*/

		REQUIRE(energy.eigenvalues     == -0.234329528903_a);
		REQUIRE(energy.xc              == -0.232294220410_a);
		REQUIRE(energy.nvxc            == -0.302713349819_a);
		REQUIRE(energy.total()         == -0.446253846698_a);
		REQUIRE(energy.external        == -0.108660738870_a);
		REQUIRE(energy.nonlocal        == -0.058633055438_a);
		REQUIRE(energy.kinetic()       ==  0.416973236003_a);
		
	}
#endif
}
