/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa

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

	using namespace inq::magnitude;

	inq::input::environment env(argc, argv);
	
	inq::utils::match energy_match(3.0e-5);

	inq::input::species local_h = pseudo::element("H") | inq::input::species::symbol("Hloc") | inq::input::species::pseudo(inq::config::path::unit_tests_data() + "H.blyp-vbc.UPF");

	auto box = inq::systems::box::cubic(15.0_b).finite().cutoff_energy(40.0_Ha);

	inq::systems::ions ions(box);

	ions.insert(local_h, {150.0_b, -30.0_b, 0.0_b});

	inq::systems::electrons electrons(env.par(), ions, box);
	inq::ground_state::initial_guess(ions, electrons);
	
	inq::ground_state::calculate(ions, electrons);
	auto result = inq::ground_state::calculate(ions, electrons, inq::input::interaction::hartree_fock(), inq::input::scf::energy_tolerance(1e-8_Ha));
	
	energy_match.check("total energy",        result.energy.total(),      -0.578525486338);
	energy_match.check("kinetic energy",      result.energy.kinetic(),     0.348185715818);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),-0.230237311450);
	energy_match.check("Hartree energy",      result.energy.hartree(),     0.254438812760);
	energy_match.check("external energy",     result.energy.external(),   -0.832861830904);
	energy_match.check("non-local energy",    result.energy.nonlocal(),    0.0);
	energy_match.check("XC energy",           result.energy.xc(),          0.0);
	energy_match.check("XC density integral", result.energy.nvxc(),        0.0);
	energy_match.check("HF exchange energy",  result.energy.hf_exchange(),-0.254438821884);
	energy_match.check("ion-ion energy",      result.energy.ion(),        -0.093849362128);
	
	return energy_match.fail();
	
}
