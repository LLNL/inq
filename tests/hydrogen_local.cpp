/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq::magnitude;

	inq::utils::match energy_match(3.0e-5);

	inq::systems::ions ions(inq::systems::cell::cubic(15.0_b).finite());
	ions.species_list().insert(inq::ionic::species("H").symbol("Hloc").pseudo_file(inq::config::path::unit_tests_data() + "H.blyp-vbc.UPF"));
	ions.insert("Hloc", {150.0_b, -30.0_b, 0.0_b});
	
	inq::systems::electrons electrons(ions, inq::options::electrons{}.cutoff(40.0_Ha));
	inq::ground_state::initial_guess(ions, electrons);
	
	// Non Interacting
	{
	
		auto result = inq::ground_state::calculate(ions, electrons, inq::options::theory{}.non_interacting(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));
		
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

		energy_match.check("total energy",        result.energy.total(),          -0.499022722967);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         0.488268502855);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.499022722967);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.000000000000);
		energy_match.check("external energy",     result.energy.external(),       -0.987291225822);
		energy_match.check("non-local energy",    result.energy.non_local(),       0.000000000000);
		energy_match.check("XC energy",           result.energy.xc(),              0.000000000000);
		energy_match.check("XC density integral", result.energy.nvxc(),            0.000000000000);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),             0.000000000000);
		
	}

	// LDA
	{
		
		auto result = inq::ground_state::calculate(ions, electrons, inq::options::theory{}.lda(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha));
		
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

		energy_match.check("total energy",        result.energy.total(),          -0.445160072256);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         0.414315464604);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.234029035766);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.281309132025);
		energy_match.check("external energy",     result.energy.external(),       -0.909266382097);
		energy_match.check("non-local energy",    result.energy.non_local(),       0.000000000000);
		energy_match.check("XC energy",           result.energy.xc(),             -0.231518286787);
		energy_match.check("XC density integral", result.energy.nvxc(),           -0.301696382322);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),             0.000000000000);
		
	}

	// B3LYP
	{

		auto result = inq::ground_state::calculate(ions, electrons, inq::options::theory{}.b3lyp(), inq::options::ground_state{}.energy_tolerance(1e-8_Ha).max_steps(300));
		
		energy_match.check("total energy",        result.energy.total(),          -0.447429725736);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         0.421657099440, 7.0e-5);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.247459996333);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.282608756325);
		energy_match.check("external energy",     result.energy.external(),       -0.916549935273, 7.0e-5);
		energy_match.check("non-local energy",    result.energy.non_local(),       0.000000000000);
		energy_match.check("XC energy",           result.energy.xc(),             -0.206884771385);
		energy_match.check("XC density integral", result.energy.nvxc(),           -0.261262923463);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(), -0.028260874843);
		energy_match.check("ion-ion energy",      result.energy.ion(),             0.000000000000);
		
	}
	
	return energy_match.fail();
	
}
