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
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
		
	utils::match energy_match(2.0e-5);

	std::vector<input::atom> geo;

	geo.push_back("Ne" | input::species::nofilter() | math::vector3<double>(0.0, 0.0, 0.0));
		
	systems::ions ions(input::cell::cubic(15.0_b) | input::cell::finite(), geo);

	input::config conf;
	
	conf.extra_states = 3;
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(30.0_Ha), conf);
	
	ground_state::initial_guess(ions, electrons);
	
	//REAL SPACE PSEUDO
	{
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting(), input::scf::conjugate_gradient());
		
		energy_match.check("total energy",     result.energy.total()    , -67.771731167338);
		energy_match.check("kinetic energy",   result.energy.kinetic()  ,  35.606535224997);
		energy_match.check("eigenvalues",      result.energy.eigenvalues, -61.765371991153);
		energy_match.check("external energy",  result.energy.external   , -79.444978433236);
		energy_match.check("non-local energy", result.energy.nonlocal   , -17.926928782914);
		energy_match.check("ion-ion energy",   result.energy.ion        , -6.006359176186);

	}

	//FOURIER SPACE PSEUDO
	{	
		auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting() | input::interaction::fourier_pseudo(), input::scf::conjugate_gradient());
		
		energy_match.check("total energy",     result.energy.total()    , -67.771735282066);
		energy_match.check("kinetic energy",   result.energy.kinetic()  ,  35.606511739929);
		energy_match.check("eigenvalues",      result.energy.eigenvalues, -61.765376105880);
		energy_match.check("external energy",  result.energy.external   , -79.444966794608);
		energy_match.check("non-local energy", result.energy.nonlocal   , -17.926921051201);
		energy_match.check("ion-ion energy",   result.energy.ion        ,  -6.006359176186);
		
	}

	return energy_match.fail();
	
}
