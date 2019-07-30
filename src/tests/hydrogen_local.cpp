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

	input::species local_h(pseudo::element("H"), config::path::unit_tests_data() + "H.blyp-vbc.UPF");
	
	ions::geometry geo;

	geo.add_atom(pseudo::element("H"), math::d3vector(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(10.0, 10.0, 10.0), geo);

	input::config conf;
	systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), conf);

	auto energy = electrons.calculate_ground_state();

	/*
	  #k =   1, k = (  0.000000,  0.000000,  0.000000)
	  1   --    -0.232923       1.000000
	  
	  Direct gap at ik=   -1 of ******* H
	  Indirect gap between ik=   -1 and ik=   -1 of     Inf H
	  
	  Energy [H]:
	  Total       =        -0.44568137
	  Free        =        -0.44568137
	  -----------
	  Ion-ion     =        -0.07062564
	  Eigenvalues =        -0.23292318
	  Hartree     =         0.21258407
	  Int[n*v_xc] =        -0.30285389
	  Exchange    =        -0.19278303
	  Correlation =        -0.03961935
	  vanderWaals =         0.00000000
	  Delta XC    =         0.00000000
	  Entropy     =         1.38629436
	  -TS         =        -0.00000000
	  Kinetic     =         0.41712467
	  External    =        -0.77236145
	  Non-local   =        -0.05816341
	*/
	
	REQUIRE(energy.total()       == -0.6848531681_a);
	REQUIRE(energy.kinetic()     == 2.368793_a);
	REQUIRE(energy.eigenvalues   == -1.6053918367_a);
	REQUIRE(fabs(energy.hartree) <=  1e-10);
	REQUIRE(energy.xc            == -3.0536456687_a);
	REQUIRE(energy.nvxc          == -3.9741843374_a);
	
}
