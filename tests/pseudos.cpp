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

	geo.push_back("Ne" | inq::input::species::pseudo(inq::config::path::unit_tests_data() + "C_ONCV_PBE-1.2.xml") | math::vector3<double>(0.0, 0.0, 0.0));
		
	systems::ions ions(systems::box::cubic(15.0_b).finite(), geo);

	input::config conf;
	
	conf.extra_states = 4;
  conf.temperature = 300.0_K;
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0_Ha), conf);
	
	ground_state::initial_guess(ions, electrons);

  auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), input::scf::conjugate_gradient());

  energy_match.check("total energy",        result.energy.total()    ,    -6.865720746337);
  energy_match.check("kinetic energy",      result.energy.kinetic()  ,     3.176990908258);
  energy_match.check("eigenvalues",         result.energy.eigenvalues,    -1.404855551387);
	energy_match.check("Hartree energy",      result.energy.hartree,         4.371246797420);
	energy_match.check("external energy",     result.energy.external,      -11.993844933515);
	energy_match.check("non-local energy",    result.energy.nonlocal,        0.494162047716);
	energy_match.check("XC energy",           result.energy.xc,             -1.412685772171);
	energy_match.check("XC density integral", result.energy.nvxc,           -1.824657168687);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,     0.000000000000);
	energy_match.check("ion-ion energy",      result.energy.ion,            -1.501589794046);
	
	return energy_match.fail();
	
}
