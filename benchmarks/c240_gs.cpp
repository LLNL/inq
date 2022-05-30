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
#include <config/path.hpp>
#include <input/atom.hpp>
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <input/parse_xyz.hpp>
#include <config/path.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;	
	
	input::environment env(argc, argv);

	auto box = systems::box::cubic(20.0_A).finite().spacing(20.0_A/90);
	
	systems::ions ions(box);

	ions.insert(input::parse_xyz(config::path::unit_tests_data() + "c240.xyz"));
	
	input::config conf;
	
	conf.extra_states = 32;
	conf.temperature = 300.0_K;
	
	systems::electrons electrons(env.par(), ions, box, conf);
	
	ground_state::initial_guess(ions, electrons);
	ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(10) | inq::input::scf::mixing(0.3));

	electrons.save("c240_restart");
}

