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
#include <input/environment.hpp>
#include <input/atom.hpp>
#include <utils/match.hpp>
#include <ground_state/calculate.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	inq::input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
		
	utils::match energy_match(1.0e-6);

	std::vector<input::atom> geo;

	geo.push_back("H" | input::species::nofilter() | math::vector3<double>(0.0, 0.0, 0.0));
		
	systems::ions ions(input::cell::cubic(10.0, 10.0, 10.0) | input::cell::finite(), geo);

	//Real space pseudo
	{
		
		input::config conf;

		conf.extra_states = 3;

		systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0_Ha), conf);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting(), input::scf::davidson());

		std::printf("total energy:%20.16f\n",     result.energy.total() );
	       
	}

	//Fourier space pseudo
	{
		
		input::config conf;

		conf.extra_states = 3;

		systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0_Ha), conf);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting() | input::interaction::fourier_pseudo(), input::scf::conjugate_gradient());

		std::printf("total energy:%20.16f\n",     result.energy.total() );
		
	}

	return 0;
	
}
