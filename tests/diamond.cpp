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
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <observables/kinetic_energy_density.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);

	utils::match energy_match(3.0e-5);

	std::vector<input::atom> geo;

	auto a =  3.567095_A;

	auto box = systems::box::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}).cutoff_energy(35.0_Ha);
	
	systems::ions ions(box);
	
	ions.insert("C", {0.0_crys,  0.0_crys,  0.0_crys });
	ions.insert("C", {0.25_crys, 0.25_crys, 0.25_crys});

	input::config conf;
	
	conf.extra_states = 3;

	systems::electrons electrons(env.par(), ions, box, conf, input::kpoints::grid({1, 1, 1}, false));

	ground_state::initial_guess(ions, electrons);

	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));

	energy_match.check("total energy",        result.energy.total(),         -10.949196617732);
	energy_match.check("kinetic energy",      result.energy.kinetic(),        11.411454690719);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,       0.795474829167);
	energy_match.check("Hartree energy",      result.energy.hartree,           1.473883973621);
	energy_match.check("external energy",     result.energy.external,         -7.034444637769);
	energy_match.check("non-local energy",    result.energy.nonlocal,         -1.496762367780);
	energy_match.check("XC energy",           result.energy.xc,               -4.568604147920);
	energy_match.check("XC density integral", result.energy.nvxc,             -5.032540803244);
	energy_match.check("ion-ion energy",      result.energy.ion,             -10.734724128603);
	
	electrons.save("diamond_restart");

	auto ked = observables::kinetic_energy_density(electrons);

	energy_match.check("kinetic energy", operations::integral(ked), 11.411455188639);
	
	fftw_cleanup();
	
	return energy_match.fail();
	
}


