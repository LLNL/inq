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

	input::species local_h = pseudo::element("H") | input::species::symbol("Hloc") | input::species::pseudo(config::path::unit_tests_data() + "H.blyp-vbc.UPF"); 
	
	std::vector<input::atom> geo;
	
	geo.push_back(local_h | math::d3vector(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);
#if 1
	SECTION("Non interacting"){
		
		input::config conf;
		
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), input::interaction::non_interacting(), conf);
		
		auto energy = electrons.calculate_ground_state();
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)
			#st  Spin   Eigenvalue      Occupation
			1   --    -0.500174       1.000000
			
			Energy [H]:
      Total       =        -0.50017433
      Free        =        -0.50017433
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =        -0.50017433
      Hartree     =         0.00000000
      Int[n*v_xc] =         0.00000000
      Exchange    =         0.00000000
      Correlation =         0.00000000
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.49296606
      External    =        -0.99314039
      Non-local   =         0.00000000

		*/
		REQUIRE(energy.self            == -0.564189583548_a);
		REQUIRE(energy.eigenvalues     == -0.499895442304_a);
		REQUIRE(fabs(energy.xc)        <=  1e-10);
		REQUIRE(fabs(energy.nvxc)      <=  1e-10);
		REQUIRE(energy.coulomb()       == -0.793703931391_a);
		REQUIRE(energy.total()         == -0.499895442304_a);
		REQUIRE(energy.external        == -0.194627483520_a);
		REQUIRE(energy.kinetic()       ==  0.488435972607_a);

	}
#endif
#if 1
	SECTION("LDA"){
		
		input::config conf;
	
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), input::interaction::dft(), conf);
		
		auto energy = electrons.calculate_ground_state();
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)

			1   --    -0.233986       1.000000

			Energy [H]:
      Total       =        -0.44606573
      Free        =        -0.44606573
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =        -0.23398591
      Hartree     =         0.28254446
      Int[n*v_xc] =        -0.30290955
      Exchange    =        -0.19282007
      Correlation =        -0.03962486
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.41903428
      External    =        -0.91520434
      Non-local   =         0.00000000

		*/

		REQUIRE(energy.self          == -0.564189583548_a); 
		REQUIRE(energy.eigenvalues   == -0.234407373371_a);
		REQUIRE(energy.xc            == -0.232064693875_a);
		REQUIRE(energy.nvxc          == -0.302412613089_a);
		REQUIRE(energy.kinetic()     ==  0.417754384376_a);
		REQUIRE(energy.total()       == -0.446068359498_a);
		REQUIRE(energy.external      == -0.168489350956_a);
		REQUIRE(energy.coulomb()     == -0.463268699041_a);
		REQUIRE(energy.ion           ==  0.000000000000_a);

	}
#endif
}
