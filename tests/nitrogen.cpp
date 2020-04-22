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

int main(int argc, char ** argv){

	boost::mpi3::environment env(argc, argv);

	utils::match energy_match(1.0e-6);

	std::vector<input::atom> geo;

	auto distance = 2.0739744;
	
	geo.push_back( "N" | math::vec3d(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | math::vec3d(0.0, 0.0,  0.5*distance));
		
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);

		input::config conf;

		conf.extra_states = 4;

		systems::electrons electrons(ions, input::basis::cutoff_energy(40.0), conf);
		
		[[maybe_unused]] auto energy = electrons.calculate_ground_state(input::interaction::dft());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)
			Eigenvalues [H]
			#st  Spin   Eigenvalue      Occupation
			1   --    -0.233880       1.000000
			
			Energy [H]:
      Total      =        -0.44618708
      Free        =        -0.44618708
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =        -0.23387966
      Hartree     =         0.28283087
      Int[n*v_xc] =        -0.30315816
      Exchange    =        -0.19299390
      Correlation =        -0.03964083
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.41796657
      External    =        -0.91434707
      Non-local   =        -0.05876129

		*/

	return energy_match.fail();
}
