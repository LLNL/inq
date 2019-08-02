/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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

TEST_CASE("Test non interacting electron gas", "[test::non_interacting_electron_gas]") {

	using namespace Catch::literals;

	input::species local_h = pseudo::element("H") | input::species::pseudo(config::path::unit_tests_data() + "H.blyp-vbc.UPF") | input::species::name("Hloc");
	
	ions::geometry geo;

	geo.add_atom(local_h, math::d3vector(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0), geo);

	//NON INTERACTING
	{
		
		input::config conf;
		
		conf.theory = input::electronic_theory::NON_INTERACTING;
		
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), conf);
		
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
		
		REQUIRE(energy.total()       == -0.4547670333_a);
		REQUIRE(energy.kinetic()     ==  0.4868395171_a);
		REQUIRE(energy.eigenvalues   == -0.3591216767_a);
		REQUIRE(energy.external      == -0.8459611938_a);
		REQUIRE(fabs(energy.hartree) <=  1e-10);
		REQUIRE(fabs(energy.xc)      <=  1e-10);
		REQUIRE(fabs(energy.nvxc)    <=  1e-10);
		REQUIRE(energy.ion           == -0.0956453566_a);
	}

	//LDA
	{
		
		input::config conf;
	
		systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), conf);
		
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
			<eself>       0.79788456 </eself>
			<ets>         0.00000000 </ets>
			<etotal>     -0.44569809 </etotal>

		*/
		
		REQUIRE(energy.total()       == -0.4712455539_a);
		REQUIRE(energy.kinetic()     ==  4.166462e-01_a);
		REQUIRE(energy.eigenvalues   == -2.342972e-01_a);
		REQUIRE(energy.external      == -7.721030e-01_a);
		REQUIRE(energy.hartree       ==  2.115318e-01_a);
		REQUIRE(energy.xc            == -2.316752e-01_a);
		REQUIRE(energy.nvxc          == -3.019039e-01_a);
		REQUIRE(energy.ion           == -0.0956453566_a);

	}
}
