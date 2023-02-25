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
	
	systems::electrons electrons(env.par(), ions, box, input::config::extra_states(1) | input::config::temperature(300.0_K), input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), input::scf::energy_tolerance(1e-8_Ha) | input::scf::calculate_forces());
	
	energy_match.check("ion-ion energy",      result.energy.ion,           -10.318372113231);
	energy_match.check("total energy",        result.energy.total(),        -9.802338589416);
	energy_match.check("kinetic energy",      result.energy.kinetic(),       4.200431340049);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,     0.602436500615);
	energy_match.check("Hartree energy",      result.energy.hartree,         0.219185817806);
	energy_match.check("external energy",     result.energy.external,       -0.562805143373);
	energy_match.check("non-local energy",    result.energy.nonlocal,        1.427216975992);
	energy_match.check("XC energy",           result.energy.xc,             -4.767995466659);
	energy_match.check("XC density integral", result.energy.nvxc,           -4.900778307666);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,     0.000000000000);
	
	energy_match.check("force 1 x",           result.forces[0][0],          -0.022483431037);
	energy_match.check("force 1 y",           result.forces[0][1],          -0.041215997171);
	energy_match.check("force 1 z",           result.forces[0][2],          -0.052723786483);
	energy_match.check("force 2 x",           result.forces[1][0],          -0.022476660700);
	energy_match.check("force 2 y",           result.forces[1][1],           0.052697035680);
	energy_match.check("force 2 z",           result.forces[1][2],           0.041207478998);
	energy_match.check("force 3 x",           result.forces[2][0],           0.005730135670);
	energy_match.check("force 3 y",           result.forces[2][1],          -0.012778476335);
	energy_match.check("force 3 z",           result.forces[2][2],           0.012775275108);
	energy_match.check("force 4 x",           result.forces[3][0],           0.007076613283);
	energy_match.check("force 4 y",           result.forces[3][1],           0.012276399154);
	energy_match.check("force 4 z",           result.forces[3][2],          -0.012280307956);
	energy_match.check("force 5 x",           result.forces[4][0],           0.027652090218);
	energy_match.check("force 5 y",           result.forces[4][1],          -0.010193515961);
	energy_match.check("force 5 z",           result.forces[4][2],           0.010356483661);
	
	electrons.save("al4h1_restart");
	
	return energy_match.fail();
	
}

