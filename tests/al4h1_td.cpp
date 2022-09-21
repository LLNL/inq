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

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);
	
	utils::match energy_match(1.0e-5);

	auto alat = 7.6524459_bohr;

	systems::box box = systems::box::cubic(alat).cutoff_energy(30.0_Ha);

	systems::ions ions(box);

	ions.insert("Al", {0.0_crys, 0.0_crys, 0.0_crys});
	ions.insert("Al", {0.0_crys, 0.5_crys, 0.5_crys});
	ions.insert("Al", {0.5_crys, 0.0_crys, 0.5_crys});
	ions.insert("Al", {0.5_crys, 0.5_crys, 0.0_crys});	
	ions.insert("H",  {0.1_crys, 0.2_crys, 0.3_crys});
	
	input::config conf;

	conf.extra_states = 1;
	conf.temperature = 300.0_K;

	systems::electrons electrons(env.par(), ions, box, conf, input::kpoints::grid({2, 2, 2}, true));
	
	electrons.load("al4h1_restart");

	std::vector<double> energy;
	
	real_time::propagate<>(ions, electrons, [&](auto data){energy.push_back(data.energy());}, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));

	energy_match.check("energy step   0", energy[0],   -9.798687545590);
	energy_match.check("energy step  10", energy[10],  -9.798687882908);
	energy_match.check("energy step  20", energy[20],  -9.798688620193);
	energy_match.check("energy step  30", energy[30],  -9.798688818773);
	
	return energy_match.fail();
	
}

