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
			
	input::environment env;
	
	utils::match energy_match(1.0e-5);

	{
		auto box = systems::cell::cubic(10.0_b);
		
		systems::ions ions(box);

        systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).temperature(300.0_K).extra_electrons(14.0).extra_states(2), input::kpoints::grid({1, 1, 1}, false));
		
		ground_state::initial_guess(ions, electrons);

        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.energy_tolerance(1e-9_Ha));
		
		energy_match.check("total energy",        result.energy.total()    , -0.567967321401);
		energy_match.check("kinetic energy",      result.energy.kinetic()  ,  2.485678165423);
//		energy_match.check("hartree",             result.energy.hartree    ,  0.000000732036);	
//		energy_match.check("XC energy",           result.energy.xc         , -3.053646218860);
	}
}
