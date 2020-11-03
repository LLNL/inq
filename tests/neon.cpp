/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019-2020 Xavier Andrade

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
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	
	boost::mpi3::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
		
	utils::match energy_match(1.0e-6);

	std::vector<input::atom> geo;

	geo.push_back("Ne" | input::species::nofilter() | math::vec3d(0.0, 0.0, 0.0));
		
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);

	//Real space pseudo
	{
		
		input::config conf;

		conf.extra_states = 3;

		systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(40.0), conf);

		ground_state::initialize(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting(), input::scf::conjugate_gradient());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.350877)

			Eigenvalues [H]
			#st  Spin   Eigenvalue      Occupation
			1   --    -8.537419       2.000000
			2   --    -7.457549       2.000000
			3   --    -7.457549       2.000000
			4   --    -7.457549       2.000000
			5   --    -3.733537       0.000000
			6   --    -3.508765       0.000000
			7   --    -3.508765       0.000000
			
			Energy [H]:
      Total       =       -61.82012973
      Free        =       -61.82012973
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =       -61.82012973
      Hartree     =         0.00000000
      Int[n*v_xc] =         0.00000000
      Exchange    =         0.00000000
      Correlation =         0.00000000
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         0.00000000
      -TS         =        -0.00000000
      Kinetic     =        35.71136464
      External    =       -97.53149437
      Non-local   =       -18.00231356

		*/
		
		energy_match.check("total energy",     result.energy.total()    , -66.260996131073);
		energy_match.check("kinetic energy",   result.energy.kinetic()  ,  35.546488765843);
		energy_match.check("eigenvalues",      result.energy.eigenvalues, -61.740955117988);
		energy_match.check("external energy",  result.energy.external   , -79.409315995073);
		energy_match.check("non-local energy", result.energy.nonlocal   , -17.878127888758);
		energy_match.check("ion-ion energy",   result.energy.ion        , -4.520041013085);
		
	}

	//Fourier space pseudo
	{
		
		input::config conf;

		conf.extra_states = 3;

		systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(40.0), conf);

		ground_state::initialize(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting() | input::interaction::fourier_pseudo(), input::scf::conjugate_gradient());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.350877)

			Eigenvalues [H]
			#st  Spin   Eigenvalue      Occupation
			1   --    -8.537419       2.000000
			2   --    -7.457549       2.000000
			3   --    -7.457549       2.000000
			4   --    -7.457549       2.000000
			5   --    -3.733537       0.000000
			6   --    -3.508765       0.000000
			7   --    -3.508765       0.000000
			
			Energy [H]:
      Total       =       -61.82012973
      Free        =       -61.82012973
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =       -61.82012973
      Hartree     =         0.00000000
      Int[n*v_xc] =         0.00000000
      Exchange    =         0.00000000
      Correlation =         0.00000000
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         0.00000000
      -TS         =        -0.00000000
      Kinetic     =        35.71136464
      External    =       -97.53149437
      Non-local   =       -18.00231356

		*/
		
		energy_match.check("total energy",     result.energy.total()    , -66.260998366598);
		energy_match.check("kinetic energy",   result.energy.kinetic()  ,  35.546473395356);
		energy_match.check("eigenvalues",      result.energy.eigenvalues, -61.740957353513);
		energy_match.check("external energy",  result.energy.external   , -79.409308401592);
		energy_match.check("non-local energy", result.energy.nonlocal   , -17.878122347276);
		energy_match.check("ion-ion energy",   result.energy.ion        ,  -4.520041013085);
		
	}

	return energy_match.fail();
	
}
