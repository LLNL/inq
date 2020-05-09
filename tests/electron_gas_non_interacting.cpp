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
#include <ground_state/calculate.hpp>

int main(int argc, char ** argv){

	boost::mpi3::environment env(argc, argv);

	utils::match energy_match(1.0e-6);

	systems::ions ions(input::cell::cubic(10.0, 10.0, 10.0));

	input::config conf;
	conf.extra_states = 2;
	conf.excess_charge = 14.0;
		
	systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), conf);

	auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting());

	//Octopus results are:
	// Energy: 2.36870506
	// Eigenvalues: 0.000000 0.197392 0.394784
	energy_match.check("total energy",   result.energy.total(),     2.3687083213);
	energy_match.check("kinetic energy", result.energy.kinetic(),   2.3687083213);
	energy_match.check("eigenvalues",    result.energy.eigenvalues, 2.3687083213);

	return energy_match.fail();

}
