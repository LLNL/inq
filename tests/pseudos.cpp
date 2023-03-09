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

	{
		systems::ions ions(box);
		
		ions.insert("C" | inq::input::species::pseudo(inq::config::path::unit_tests_data() + "C_ONCV_PBE-1.2.xml"), {0.0_b, 0.0_b, 0.0_b});

		systems::electrons electrons(env.par(), ions, box, input::config::extra_states(4) | input::config::temperature(300.0_K));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total()    ,   -5.364150140372);
		energy_match.check("kinetic energy",      result.energy.kinetic()   ,   3.176958479664);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(), -1.404891046908);
		energy_match.check("Hartree energy",      result.energy.hartree(),      4.371229020227);
		energy_match.check("external energy",     result.energy.external(),   -11.993822927726);
		energy_match.check("non-local energy",    result.energy.nonlocal(),     0.494166478064);
		energy_match.check("XC energy",           result.energy.xc(),          -1.412681190601);
		energy_match.check("XC density integral", result.energy.nvxc(),        -1.824651117364);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),          0.000000000000);

	}

	{
		systems::ions ions(box);
		
		ions.insert("C" | inq::input::species::pseudo(inq::config::path::unit_tests_data() + "C.ccECP.upf"), {0.0_b, 0.0_b, 0.0_b});

		systems::electrons electrons(env.par(), ions, box, input::config::extra_states(4) | input::config::temperature(300.0_K));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total()    ,   -5.390943623548);
		energy_match.check("kinetic energy",      result.energy.kinetic()  ,    3.462168739143);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(), -1.409268764339);
		energy_match.check("Hartree energy",      result.energy.hartree(),      4.396047239823);
		energy_match.check("external energy",     result.energy.external(),   -12.434345771131);
		energy_match.check("non-local energy",    result.energy.nonlocal(),     0.607257611969);
		energy_match.check("XC energy",           result.energy.xc(),          -1.422071443352);
		energy_match.check("XC density integral", result.energy.nvxc(),        -1.836443823966);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),  0.000000000000);
		energy_match.check("ion-ion energy",      result.energy.ion(),          0.000000000000);
	}
	return energy_match.fail();
	
}
