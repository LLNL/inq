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
	
	utils::match energy_match(1.0e-6);

	auto box = systems::box::cubic(10.0_b).cutoff_energy(40.0_Ha);
	
	systems::ions ions(box);

	input::config conf;
	conf.extra_states = 5;
	conf.excess_charge = 14.0;

	systems::electrons electrons(env.par(), ions, box, conf);

	ground_state::initial_guess(ions, electrons);
	auto result = ground_state::calculate(ions, electrons, input::interaction::lda());

	energy_match.check("total energy",        result.energy.total()    , -0.684795414294);
	energy_match.check("kinetic energy",      result.energy.kinetic()  ,  2.368849520758);
	energy_match.check("eigenvalues",         result.energy.eigenvalues, -1.605333033837);
	energy_match.check("eigenvalues",         result.energy.hartree    ,  0.000001191289);	
	energy_match.check("XC energy",           result.energy.xc         , -3.053646126342);
	energy_match.check("XC density integral", result.energy.nvxc       , -3.974184937174);

	return energy_match.fail();
	
}
