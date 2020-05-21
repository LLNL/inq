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

#include <fftw3.h>

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

	utils::match energy_match(2.5e-4);

	std::vector<input::atom> geo;

	geo.push_back( "O" | math::vec3d( 0.0,      -0.553586, 0.0));
	geo.push_back( "H" | math::vec3d( 1.429937,  0.553586, 0.0));
  geo.push_back( "H" | math::vec3d(-1.429937,  0.553586, 0.0));

	systems::ions ions(input::cell::cubic(12.0, 11.0, 10.0) | input::cell::finite(), geo);

  input::config conf;
  
  systems::electrons electrons(ions, input::basis::cutoff_energy(20.0), conf);

	operations::io::load("h2o_restart", electrons.phi_);

  perturbations::kick({1.0, 0.0, 0.0}, electrons.phi_);
											 
  auto result = real_time::propagate(electrons);

  energy_match.check("energy step  0", result.energy[0],  -16.678775673923);
  energy_match.check("energy step 10", result.energy[10], -16.679224520699);

  {
    auto dipole_file = std::ofstream("dipole.dat");

    for(unsigned ii = 0; ii < result.dipole.size(); ii++){
      dipole_file << ii << '\t' << result.time[ii] << '\t' << result.dipole[ii] << std::endl;
    }
  }

	fftw_cleanup(); //required for valgrid
	
	return energy_match.fail();
}
