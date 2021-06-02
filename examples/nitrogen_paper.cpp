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
	auto comm = boost::mpi3::environment::get_world_instance();

	auto distance = 1.06_angstrom;

	std::vector<atom> geo;
	geo.push_back( "N" | coord(0.0, 0.0, -distance/2));
	geo.push_back( "N" | coord(0.0, 0.0,  distance/2));
	
	systems::ions ions(cell::orthorhombic(10.0_b, 10.0_b, 12.0_b) | cell::finite(), geo);

	systems::electrons electrons(comm, ions, input::basis::cutoff_energy(40.0_Ha));
	ground_state::initialize(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, interaction::pbe());

	std::cout << "N2 energy = " << result.energy.total() << std::endl;

}

