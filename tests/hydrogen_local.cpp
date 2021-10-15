/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa

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
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq::magnitude;

	inq::input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	inq::utils::match energy_match(3.0e-5);

	inq::input::species local_h = pseudo::element("H") | inq::input::species::symbol("Hloc") | inq::input::species::pseudo(inq::config::path::unit_tests_data() + "H.blyp-vbc.UPF");
	
	std::vector<inq::input::atom> geo;
	
	geo.push_back(local_h | inq::math::vector3<double>(150.0, -30.0, 0.0));

	inq::systems::ions ions(inq::systems::box::cubic(15.0_b) | inq::systems::box::finite(), geo);

	inq::input::config conf;
	
	inq::systems::electrons electrons(comm_world, ions, inq::input::basis::cutoff_energy(40.0_Ha), conf);
	inq::ground_state::initial_guess(ions, electrons);
	
	// Non Interacting
	{
	
		auto result = inq::ground_state::calculate(ions, electrons, inq::input::interaction::non_interacting(), inq::input::scf::conjugate_gradient() | inq::input::scf::energy_tolerance(1e-8_Ha));
		
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

		energy_match.check("total energy",        result.energy.total(),      -0.592872083038);
		energy_match.check("kinetic energy",      result.energy.kinetic(),     0.488259135290);
		energy_match.check("eigenvalues",         result.energy.eigenvalues,  -0.499022720910);
		energy_match.check("Hartree energy",      result.energy.hartree,       0.000000000000);
		energy_match.check("external energy",     result.energy.external,     -0.987281856200);
		energy_match.check("non-local energy",    result.energy.nonlocal,      0.0);
		energy_match.check("XC energy",           result.energy.xc,            0.0);
		energy_match.check("XC density integral", result.energy.nvxc,          0.0);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange,   0.0);
		energy_match.check("ion-ion energy",      result.energy.ion,          -0.093849362128);
		
	}

	// LDA
	{
		
		auto result = inq::ground_state::calculate(ions, electrons);
		
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

		energy_match.check("total energy",        result.energy.total(),          -0.539007053491);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         0.414318995673);
		energy_match.check("eigenvalues",         result.energy.eigenvalues,      -0.234026421108);
		energy_match.check("Hartree energy",      result.energy.hartree,           0.281309443919);
		energy_match.check("external energy",     result.energy.external,         -0.909267595449);
		energy_match.check("non-local energy",    result.energy.nonlocal,          0.0);
		energy_match.check("XC energy",           result.energy.xc,               -0.231518535506);
		energy_match.check("XC density integral", result.energy.nvxc,             -0.301696709170);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange,       0.0);
		energy_match.check("ion-ion energy",      result.energy.ion,              -0.093849362128);
		
	}
	
	return energy_match.fail();
	
}
