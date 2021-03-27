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

	using namespace inq;
	using namespace inq::magnitude;	
	
	input::environment env(argc, argv);

	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(4.0e-6);

	auto geo = input::parse_xyz(config::path::unit_tests_data() + "al256.xyz");

	systems::ions ions(input::cell::cubic(4*7.6524459_b), geo);
	
	input::config conf;
	
	conf.extra_states = 0;
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0_Ha), conf);
	
	ground_state::initialize(ions, electrons);

	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(2));
	
	energy_match.check("total energy",        result.energy.total(),         -1920.807756665174);
	energy_match.check("kinetic energy",      result.energy.kinetic(),         327.618032246892);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,       424.360184739629);
	energy_match.check("Hartree energy",      result.energy.hartree,          1580.487552779374);
	energy_match.check("external energy",     result.energy.external,        -2730.886670303735);
	energy_match.check("non-local energy",    result.energy.nonlocal,           49.326855481987);
	energy_match.check("XC energy",           result.energy.xc,               -356.597956355480);
	energy_match.check("XC density integral", result.energy.nvxc,             -382.673138244262);
	energy_match.check("ion-ion energy",      result.energy.ion,              -790.755570514212);
	
	fftw_cleanup(); //required for valgrid
	
	return energy_match.fail();
	
}

