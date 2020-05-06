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

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <utils/match.hpp>
#include <ground_state/calculate.hpp>

#include <mpi3/environment.hpp>

int main(int argc, char ** argv){

	boost::mpi3::environment env(argc, argv);

	utils::match energy_match(1.0e-5);

	input::species local_h = pseudo::element("H") | input::species::symbol("Hloc") | input::species::pseudo(config::path::unit_tests_data() + "H.blyp-vbc.UPF"); 
	
	std::vector<input::atom> geo;
	
	geo.push_back(local_h | math::vec3d(0.0, 0.0, 0.0));
    
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);

	input::config conf;
	
	systems::electrons electrons(ions, input::basis::cutoff_energy(60.0), conf);

	// Non Interacting
	{
	
		auto energy = ground_state::calculate(electrons, input::interaction::non_interacting(), input::scf::conjugate_gradient());
		
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

		energy_match.check("total energy",        energy.total(),      -0.570687439533);
		energy_match.check("kinetic energy",      energy.kinetic(),     0.492500951533);
		energy_match.check("eigenvalues",         energy.eigenvalues,  -0.500061798704);
		energy_match.check("Hartree energy",      energy.hartree,       0.0);		
		energy_match.check("external energy",     energy.external,     -0.992562750237);
		energy_match.check("non-local energy",    energy.nonlocal,      0.0);
		energy_match.check("XC energy",           energy.xc,            0.0);
		energy_match.check("XC density integral", energy.nvxc,          0.0);
		energy_match.check("HF exchange energy",  energy.hf_exchange,   0.0);
		energy_match.check("ion-ion energy",      energy.ion,          -0.070625640829);
		
	}

	// LDA
	{
		
		auto energy = ground_state::calculate(electrons);
		
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

		energy_match.check("total energy",        energy.total(),          -0.516602958211);
		//octopus                                                           0.41903428
		energy_match.check("kinetic energy",      energy.kinetic(),         0.418581442489);
		//octopus                                                          -0.23398591
		energy_match.check("eigenvalues",         energy.eigenvalues,      -0.233981484327);
		//octopus                                                           0.28254446
		energy_match.check("Hartree energy",      energy.hartree,           0.282434181127);
		//octopus                                                          -0.91520434
		energy_match.check("external energy",     energy.external,         -0.914632814228);
		energy_match.check("non-local energy",    energy.nonlocal,          0.0);
		//octopus                                                          -0.23244493
		energy_match.check("XC energy",           energy.xc,               -0.232360126770);
		//octopus                                                          -0.30290955
		energy_match.check("XC density integral", energy.nvxc,             -0.302798474842);
		energy_match.check("HF exchange energy",  energy.hf_exchange,       0.0);
		energy_match.check("ion-ion energy",      energy.ion,              -0.070625640829);
		
	}

	return energy_match.fail();
	
}
