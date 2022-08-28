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
		
	utils::match energy_match(2.0e-5);

	auto box = systems::box::cubic(15.0_b).finite().cutoff_energy(25.0_Ha);
	
	systems::ions ions(box);

	ions.insert("Ne" | inq::input::species::pseudo(inq::config::path::unit_tests_data() + "C_ONCV_PBE-1.2.xml"), {0.0_b, 0.0_b, 0.0_b});

	input::config conf;
	
	conf.extra_states = 4;
  conf.temperature = 300.0_K;
	
	systems::electrons electrons(env.par(), ions, box, conf);
	
	ground_state::initial_guess(ions, electrons);

  auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::energy_tolerance(1e-7_Ha));

  energy_match.check("total energy",        result.energy.total()    ,   -5.364150410132);
  energy_match.check("kinetic energy",      result.energy.kinetic()  ,    3.176958399971);
  energy_match.check("eigenvalues",         result.energy.eigenvalues,   -1.404891089837);
	energy_match.check("Hartree energy",      result.energy.hartree,        4.371228814208);
	energy_match.check("external energy",     result.energy.external,      -11.993822537719);
	energy_match.check("non-local energy",    result.energy.nonlocal,        0.494166459143);
	energy_match.check("XC energy",           result.energy.xc,             -1.412681130237);
	energy_match.check("XC density integral", result.energy.nvxc,           -1.824651039648);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange,     0.000000000000);
	energy_match.check("ion-ion energy",      result.energy.ion,             0.000000000000);
	
	return energy_match.fail();
	
}
