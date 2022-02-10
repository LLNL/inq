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
	
	ions.insert("Si", a*math::vector3<double>{0.0,  0.0,  0.0 });
	ions.insert("Si", a*math::vector3<double>{0.25, 0.25, 0.25});
	ions.insert("Si", a*math::vector3<double>{0.5,  0.5,  0.0 });
	ions.insert("Si", a*math::vector3<double>{0.75, 0.75, 0.25});
	ions.insert("Si", a*math::vector3<double>{0.5,  0.0,  0.5 });
	ions.insert("Si", a*math::vector3<double>{0.75, 0.25, 0.75});
	ions.insert("Si", a*math::vector3<double>{0.0,  0.5,  0.5 });
	ions.insert("Si", a*math::vector3<double>{0.25, 0.75, 0.75});
	
	input::config conf;
	
	conf.extra_states = 4;
	
	systems::electrons electrons(env.par(), ions, box, conf);

	auto result = real_time::propagate<>(ions, electrons, input::interaction::non_interacting(), input::rt::num_steps(100) | input::rt::dt(0.055_atomictime));
	
	electrons.load("silicon_restart");

	return energy_match.fail();
	
}

