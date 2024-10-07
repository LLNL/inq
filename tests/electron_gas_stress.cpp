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

template<class LatticeParameters>
auto total_energy(inq::input::environment& env, LatticeParameters lat) {
	using namespace inq;
	using namespace inq::magnitude;

    auto const box = systems::cell::cubic(lat);

    systems::ions const ions(box);

    systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).temperature(300.0_K).extra_electrons(14.0).extra_states(2), input::kpoints::grid({1, 1, 1}, false));

    ground_state::initial_guess(ions, electrons);

    return ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.energy_tolerance(1e-9_Ha));
}

int main(int argc, char** argv){

	using namespace inq;
	using namespace inq::magnitude;
			
	input::environment env;
	
	utils::match energy_match(1.0e-5);

    auto const lat1_in_bohr = 10.0;
    auto const vol1_in_bohr3 = std::pow(lat1_in_bohr, 3);

    auto const result1 = total_energy(env, lat1_in_bohr*1.0_b);

    energy_match.check("total energy",        result1.energy.total()    , -0.567967321401);
    energy_match.check("kinetic energy",      result1.energy.kinetic()  ,  2.485678165423);

    auto const lat2_in_bohr = 10.01;
    auto const vol2_in_bohr3 = std::pow(lat2_in_bohr, 3);

    auto const result2 = total_energy(env, lat2_in_bohr*1.0_b);
}
