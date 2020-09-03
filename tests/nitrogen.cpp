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

int main(int argc, char ** argv){

	using namespace inq;
	
	boost::mpi3::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(1.0e-6);

	std::vector<input::atom> geo;

	auto distance = 2.0739744;
	
	geo.push_back( "N" | math::vec3d(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | math::vec3d(0.0, 0.0,  0.5*distance));
		
	systems::ions ions(input::cell::cubic(20.0, 20.0, 20.0) | input::cell::finite(), geo);

		input::config conf;

		conf.extra_states = 4;

		systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(40.0), conf);

		ground_state::initialize(ions, electrons);
		[[maybe_unused]] auto result = ground_state::calculate(ions, electrons, input::interaction::dft(), input::scf::conjugate_gradient() | input::scf::density_mixing());
		
		/*
			OCTOPUS RESULTS: (Spacing 0.286)

		*/
		
		energy_match.check("ion-ion energy",      result.energy.ion,               5.020189258245);
		energy_match.check("eigenvalues",         result.energy.eigenvalues,      -5.595065831555);
		energy_match.check("total energy",        result.energy.total(),         -27.666250041632);
		energy_match.check("kinetic energy",      result.energy.kinetic(),        13.458308963383);
		energy_match.check("Hartree energy",      result.energy.hartree,          27.748075671835);
		energy_match.check("external energy",     result.energy.external,        -66.394506174992);
		energy_match.check("non-local energy",    result.energy.nonlocal,         -1.706461598696);
		energy_match.check("XC energy",           result.energy.xc,               -5.791856161407);
		energy_match.check("XC density integral", result.energy.nvxc,             -6.448558364920);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange,       0.0);

		
	return energy_match.fail();
}
