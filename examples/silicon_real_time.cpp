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
#include <real_time/propagate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);
	
	utils::match energy_match(1.0e-6);

	auto a = 10.18_bohr;

	auto box = systems::box::cubic(a).cutoff_energy(25.0_Ha);
	systems::ions ions(box);
	
	ions.insert("Si", {0.0_crys,  0.0_crys,  0.0_crys });
	ions.insert("Si", {0.25_crys, 0.25_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.5_crys,  0.0_crys });
	ions.insert("Si", {0.75_crys, 0.75_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.0_crys,  0.5_crys });
	ions.insert("Si", {0.75_crys, 0.25_crys, 0.75_crys});
	ions.insert("Si", {0.0_crys,  0.5_crys,  0.5_crys });
	ions.insert("Si", {0.25_crys, 0.75_crys, 0.75_crys});
	
	input::config conf;
	
	conf.extra_states = 4;
	
	systems::electrons electrons(env.par(), ions, box, conf);

	real_time::propagate<>(ions, electrons, [](auto){}, input::interaction::non_interacting(), input::rt::num_steps(100) | input::rt::dt(0.055_atomictime));
	
	electrons.load("silicon_restart");

	return energy_match.fail();
	
}

