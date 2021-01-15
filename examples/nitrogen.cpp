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
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	
	input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(1.0e-6);

	std::vector<input::atom> geo;

	auto distance = 2.2; //a bit larger than experiment to check the force
	
	geo.push_back( "N" | math::vector3<double>(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | math::vector3<double>(0.0, 0.0,  0.5*distance));
		
	systems::ions ions(input::cell::cubic(10.0, 10.0, 12.0), geo);

	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0), input::config{});
	ground_state::initialize(ions, electrons);
	
	auto result = ground_state::calculate(ions, electrons, input::interaction::dft(), input::scf::steepest_descent() | input::scf::density_mixing());
	
	energy_match.check("ion-ion energy",      result.energy.ion,             -1.517434464849);

	energy_match.check("total energy",        result.energy.total(),         -20.642638450082);
	energy_match.check("kinetic energy",      result.energy.kinetic(),        13.163479187955);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,      -5.266545183839);
	energy_match.check("Hartree energy",      result.energy.hartree,          14.494239370982);
	energy_match.check("external energy",     result.energy.external,        -39.444396029142);
	energy_match.check("non-local energy",    result.energy.nonlocal,         -1.620128453793);
	energy_match.check("XC energy",           result.energy.xc,               -5.718398061234);
	energy_match.check("XC density integral", result.energy.nvxc,             -6.353978630822);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,       0.0);
	
	energy_match.check("force 1 x",           result.forces[0][0],            -0.000000025966);
	energy_match.check("force 1 y",           result.forces[0][1],             0.000000017231);
	energy_match.check("force 1 z",           result.forces[0][2],             0.145730384341);
	energy_match.check("force 2 x",           result.forces[1][0],             0.000000035122);
	energy_match.check("force 2 y",           result.forces[1][1],            -0.000000025695);
	energy_match.check("force 2 z",           result.forces[1][2],            -0.145731033675);
	
	auto dt = 0.055;
	
	auto propagation = real_time::propagate(ions, electrons, input::interaction::non_interacting(), input::rt::num_steps(1000) | input::rt::dt(dt), ions::propagator::impulsive{});

	return energy_match.fail();

}

