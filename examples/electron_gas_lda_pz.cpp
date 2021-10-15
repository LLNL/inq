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

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
			
	input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(1.0e-6);

	auto box = systems::box::cubic(10.0_b).cutoff_energy(40.0_Ha);
	
	systems::ions ions(box);

	input::config conf;
	conf.extra_states = 2;
	conf.excess_charge = 14.0;

	systems::electrons electrons(comm_world, ions, box, conf);

	ground_state::initial_guess(ions, electrons);
	auto result = ground_state::calculate(ions, electrons, input::interaction::dft());

	/* OCTOPUS RESULTS (Spacing = 0.2):

		 Eigenvalues [H]
		 1   --    -0.283870       2.000000
		 2   --    -0.086478       2.000000
		 3   --    -0.086478       2.000000
		 4   --    -0.086478       2.000000
		 5   --    -0.086478       2.000000
		 6   --    -0.086478       2.000000
		 7   --    -0.086478       2.000000
		 8   --     0.110914       0.000000
		 9   --     0.110914       0.000000

		 Energy [H]:
		 Total       =       -28.66787951
		 Free        =       -28.66787951
		 -----------
		 Ion-ion     =       -27.98293976
		 Eigenvalues =        -1.60547842
		 Hartree     =         0.00000000
		 Int[n*v_xc] =        -3.97418434
		 Exchange    =        -2.49204438
		 Correlation =        -0.56160129
		 vanderWaals =         0.00000000
		 Delta XC    =         0.00000000
		 Entropy     =         0.00000000
		 -TS         =        -0.00000000
		 Kinetic     =         2.36870515
		 External    =         0.00000000
		 Non-local   =         0.00000000
	*/

	energy_match.check("total energy",        result.energy.total()    , -0.6848531681);
	energy_match.check("kinetic energy",      result.energy.kinetic()  , 2.368793);
	energy_match.check("eigenvalues",         result.energy.eigenvalues, -1.6053918367);
	energy_match.check("XC energy",           result.energy.xc         , -3.0536456687);
	energy_match.check("XC density integral", result.energy.nvxc       , -3.9741843374);

	return energy_match.fail();
	
}
