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

TEST_CASE("Test non interacting electron gas", "[test::non_interacting_electron_gas]") {

	using namespace Catch::literals;

	input::species local_h = pseudo::element("H") | input::species::symbol("Hloc") | input::species::pseudo(config::path::unit_tests_data() + "H.blyp-vbc.UPF"); 
	
	std::vector<input::atom> geo;
	
	geo.push_back(local_h | math::d3vector(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0), geo);

	SECTION("Non interacting"){
		
		input::config conf;
		
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), input::interaction::non_interacting(), conf);
		
		auto energy = electrons.calculate_ground_state();
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)

			   1   --    -0.359408       1.000000
				 
				 Total       =        -0.43003364
				 Free        =        -0.43003364
				 -----------
				 Ion-ion     =        -0.07062564
				 Eigenvalues =        -0.35940800
				 Hartree     =         0.00000000
				 Int[n*v_xc] =         0.00000000
				 Exchange    =         0.00000000
				 Correlation =         0.00000000
				 vanderWaals =         0.00000000
				 Delta XC    =         0.00000000
				 Entropy     =         1.38629436
				 -TS         =        -0.00000000
				 Kinetic     =         0.49139567
				 External    =        -0.85080367
				 Non-local   =         0.00000000

		*/
		REQUIRE(energy.self            == -0.7978845608_a); //this matches qball exactly
		REQUIRE(energy.eigenvalues     == -0.3590112743_a);
		REQUIRE(fabs(energy.xc)        <=  1e-10);
		REQUIRE(fabs(energy.nvxc)      <=  1e-10);
		REQUIRE(fabs(energy.coulomb()) <=  1e-10);
		REQUIRE(energy.total()         == -0.4547670333_a);
		REQUIRE(energy.external        == -0.8459611938_a);
		REQUIRE(energy.ion             == -0.0956453566_a);
	}
#if 1
	SECTION("LDA"){
		
		input::config conf;
	
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), input::interaction::dft(), conf);
		
		auto energy = electrons.calculate_ground_state();
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)

			1   --    -0.233333       1.000000

      Total       =        -0.44607691
      Free        =        -0.44607691
      -----------
      Ion-ion     =        -0.07062564
      Eigenvalues =        -0.23333336
      Hartree     =         0.21257486
      Int[n*v_xc] =        -0.30287671
      Exchange    =        -0.19279832
      Correlation =        -0.03962144
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.41896742
      External    =        -0.77457392
      Non-local   =         0.00000000

			
			QBALL RESULTS

			<eigenvalue_sum> -0.23303692 </eigenvalue_sum>

			<ekin>        0.41725484 </ekin>
			<econf>       0.00000000 </econf>
			<eps>        -0.10437132 </eps>
			<enl>         0.00000000 </enl>
			<ecoul>      -0.52645873 </ecoul>
			<exc>        -0.23212287 </exc>
			<evdw>        0.00000000 </evdw>
			<esr>         0.00000000 </esr>
			<eself>       0.79788456 </eself> //this is -self.
			<ets>         0.00000000 </ets>
			<etotal>     -0.44569809 </etotal>

		*/

		REQUIRE(energy.self          == -0.7978845608_a); //this matches qball exactly
		REQUIRE(energy.eigenvalues   == -0.2334309275_a); //this is a reasonable match, but I would expect it to be closer to qball since the calculation is essentially the same
		REQUIRE(energy.xc            == -0.2321787884_a);
		REQUIRE(energy.nvxc          == -0.3025616239_a);
		REQUIRE(energy.total()       == -0.4712455539_a);
		REQUIRE(energy.external      == -7.721030e-01_a);
		REQUIRE(energy.coulomb()     ==  2.115318e-01_a);
		REQUIRE(energy.ion           == -0.0956453566_a);

	}
#endif
}
