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

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {boost::mpi3::fill, 1});
	
	utils::match energy_match(4.0e-6);

	double alat = 7.6524459;
	
	std::vector<input::atom> cell;

	cell.push_back( "Al" | alat*math::vector3<double>(0.0, 0.0, 0.0));
	cell.push_back( "Al" | alat*math::vector3<double>(0.0, 0.5, 0.5));
	cell.push_back( "Al" | alat*math::vector3<double>(0.5, 0.0, 0.5));
	cell.push_back( "Al" | alat*math::vector3<double>(0.5, 0.5, 0.0));	
	cell.push_back( "H"  | alat*math::vector3<double>(0.1, 0.2, 0.3));
	
	systems::box box = systems::box::orthorhombic(alat*1.0_b, alat*1.0_b, alat*1.0_b).cutoff_energy(30.0_Ha);
	systems::ions ions(box, cell);
	
	input::config conf;

	conf.extra_states = 1;
	conf.temperature = 300.0_K;	
	
	systems::electrons electrons(cart_comm, ions, box, conf, input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe());
	
	return energy_match.fail();
	
}

