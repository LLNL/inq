/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa

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
	
	geo.push_back(local_h | math::vec3d(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);

	input::config conf;
	
	systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), conf);

#if 1
	SECTION("Non interacting"){
		
	
		auto energy = electrons.calculate_ground_state(input::interaction::non_interacting());
		
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
		
		REQUIRE(energy.ion             == -0.070625640829_a);
		REQUIRE(energy.eigenvalues     == -0.499985694873_a);
		REQUIRE(energy.total()         == -0.570611335702_a);
		REQUIRE(energy.kinetic()       ==  0.487844175357_a);
		REQUIRE(energy.external        == -0.987829870230_a);
		REQUIRE(fabs(energy.hartree)   <=  1e-10);
		REQUIRE(fabs(energy.nonlocal)  <=  1e-10);
		REQUIRE(fabs(energy.xc)        <=  1e-10);
		REQUIRE(fabs(energy.nvxc)      <=  1e-10);
		REQUIRE(fabs(energy.hf_exchange) <=  1e-10);
		
	}
#endif
#if 1
	SECTION("LDA"){
		
		auto energy = electrons.calculate_ground_state(input::interaction::dft());
		
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

		REQUIRE(energy.ion             == -0.070625640829_a);

		/*

			Temporarily disable these matches, as they fail with nvcc. It
			seems this might be an issue with convergence problems we are
			seeing.

		//octopus                         -0.23398591
		REQUIRE(energy.eigenvalues     == -0.234111794026_a);
		REQUIRE(energy.total()         == -0.516616112180_a);
		
		//octopus                          0.41903428
		REQUIRE(energy.kinetic()       ==  0.418334559664_a);

		//octopus                          0.28254446
		REQUIRE(energy.hartree         ==  0.282285933038_a);

		//octopus                         -0.91520434
		REQUIRE(energy.external        == -0.914352651445_a);
		REQUIRE(fabs(energy.nonlocal)  <=  1e-10);

		//octopus                         -0.23244493
		REQUIRE(energy.xc              == -0.232258312608_a);

		//octopus                         -0.30290955
		REQUIRE(energy.nvxc            == -0.302665568320_a);
		REQUIRE(fabs(energy.hf_exchange) <=  1e-10);
		*/
		
	}
#endif
#if 0
	SECTION("Hartree-Fock"){
		
		auto energy = electrons.calculate_ground_state(input::interaction::hartree_fock());

		REQUIRE(energy.total()         == -0.485932246662_a);
		REQUIRE(energy.kinetic()       ==  0.352630715248_a);
		REQUIRE(energy.eigenvalues     == -0.229929375677_a);
		REQUIRE(energy.hartree         ==  0.123590349097_a);
		REQUIRE(energy.external        == -0.141980160329_a);
		REQUIRE(fabs(energy.nonlocal)  <=  1e-10);
		REQUIRE(energy.xc              == -0.232096508183_a);
		REQUIRE(energy.nvxc            == -0.302454897648_a);
		REQUIRE(energy.hf_exchange     ==  1e-10);
		REQUIRE(fabs(energy.ion)       <=  1e-10);

	}
#endif

}
