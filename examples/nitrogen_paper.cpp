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
	using namespace inq::input;
	using namespace inq::magnitude;
	
	environment env(argc, argv);

	auto distance = 1.06_angstrom;

	auto box = systems::box::orthorhombic(10.0_b, 10.0_b, 12.0_b).finite().cutoff_energy(40.0_Ha);
	
	systems::ions ions(box);

	ions.insert("N", {0.0_b, 0.0_b, -distance/2});
  ions.insert("N", {0.0_b, 0.0_b,  distance/2});
	
	systems::electrons electrons(env.par(), ions, box);
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, interaction::pbe());

	std::cout << "N2 energy = " << result.energy.total() << std::endl;

}

