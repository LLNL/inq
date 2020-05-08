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
#include <operations/io.hpp>
#include <perturbations/kick.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>

int main(int argc, char ** argv){

	boost::mpi3::environment env(argc, argv);

	utils::match match(1e-5);

	std::vector<input::atom> geo;

	geo.push_back( "O" | math::vec3d( 0.0,      -0.553586, 0.0));
	geo.push_back( "H" | math::vec3d( 1.429937,  0.553586, 0.0));
  geo.push_back( "H" | math::vec3d(-1.429937,  0.553586, 0.0));

	systems::ions ions(input::cell::cubic(12.0, 11.0, 10.0) | input::cell::finite(), geo);

  input::config conf;
  
  systems::electrons electrons(ions, input::basis::cutoff_energy(20.0), conf);

	// Propagation without perturbation
	{
		operations::io::load("h2o_restart", electrons.phi_);
		
		auto result = real_time::propagate(electrons, input::interaction::dft(), input::rt::num_steps(100) | input::rt::dt(0.055));
		
		match.check("energy step   0", result.energy[0],   -25.885010471387);
		match.check("energy step  10", result.energy[10],  -25.885010471211);
		match.check("energy step  20", result.energy[20],  -25.885010471034);
		match.check("energy step  30", result.energy[30],  -25.885010470857);
		match.check("energy step  40", result.energy[40],  -25.885010470679);
		match.check("energy step  50", result.energy[50],  -25.885010470501);
		match.check("energy step  60", result.energy[60],  -25.885010470322);
		match.check("energy step  70", result.energy[70],  -25.885010470143);
		match.check("energy step  80", result.energy[80],  -25.885010469964);
		match.check("energy step  90", result.energy[90],  -25.885010469785);
		match.check("energy step 100", result.energy[100], -25.885010469605);

	}

	/*
		
	{
		operations::io::load("h2o_restart", electrons.phi_);
		
		perturbations::kick({0.1, 0.0, 0.0}, electrons.phi_);
		
		auto result = real_time::propagate(electrons);
		
		match.check("energy step  0", result.energy[0], -16.903925978590);
		match.check("energy step 10", result.energy[10], -16.904635586794);
		
		{
			auto dipole_file = std::ofstream("dipole.dat");
			
			for(unsigned ii = 0; ii < result.dipole.size(); ii++){
				dipole_file << ii << '\t' << result.time[ii] << '\t' << result.dipole[ii] << std::endl;
			}
		}
	}
	*/
	return match.fail();
}
