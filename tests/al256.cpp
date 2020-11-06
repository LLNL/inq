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
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <input/parse_xyz.hpp>
#include <config/path.hpp>

int main(int argc, char ** argv){

	CALI_CXX_MARK_FUNCTION;

	using namespace inq;
	
	input::environment env(argc, argv);

	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(4.0e-6);

	auto geo = input::parse_xyz(config::path::unit_tests_data() + "al256.xyz");

	systems::ions ions(input::cell::cubic(4*7.6524459), geo);
	
	input::config conf;
	
	conf.extra_states = 0;
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0), conf);
	
	ground_state::initialize(ions, electrons);

	//auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(2));
	
	/*	energy_match.check("total energy",     result.energy.total()    , -23.695217119917);
	energy_match.check("kinetic energy",   result.energy.kinetic()  ,  14.889578471466);
	energy_match.check("eigenvalues",      result.energy.eigenvalues,   7.788403571539);
	energy_match.check("external energy",  result.energy.external   , -12.295872226537);
	energy_match.check("non-local energy", result.energy.nonlocal   ,   5.194697326610);
	energy_match.check("ion-ion energy",   result.energy.ion        , -31.483620691456);
	*/

	fftw_cleanup(); //required for valgrid
	
	return energy_match.fail();
	
}

