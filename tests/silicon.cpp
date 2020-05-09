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
#include <utils/match.hpp>
#include <ground_state/calculate.hpp>

int main(int argc, char ** argv){

	boost::mpi3::environment env(argc, argv);

	utils::match energy_match(1.0e-6);

	std::vector<input::atom> geo;

	double a = 10.18;
	
	geo.push_back( "Si" | a*math::vec3d(0.0,  0.0,  0.0 ));
	geo.push_back( "Si" | a*math::vec3d(0.25, 0.25, 0.25));
	geo.push_back( "Si" | a*math::vec3d(0.5,  0.5,  0.0 ));
	geo.push_back( "Si" | a*math::vec3d(0.75, 0.75, 0.25));
	geo.push_back( "Si" | a*math::vec3d(0.5,  0.0,  0.5 ));
	geo.push_back( "Si" | a*math::vec3d(0.75, 0.25, 0.75));
	geo.push_back( "Si" | a*math::vec3d(0.0,  0.5,  0.5 ));
	geo.push_back( "Si" | a*math::vec3d(0.25, 0.75, 0.75));

	systems::ions ions(input::cell::cubic(a) | input::cell::finite(), geo);
	
	input::config conf;
	
	conf.extra_states = 4;
	
	systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), conf);
	
	[[maybe_unused]] auto result = ground_state::calculate(ions, electrons, input::interaction::dft());
	
	/*
		OCTOPUS RESULTS: (Spacing 0.286)
	*/
		
}

