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

	auto dcc = 1.42_A;
  auto aa = sqrt(3)*dcc;
  auto lz = 10.0_b;

	auto box = systems::box::lattice(aa*vector3{1.0, 0.0, 0.0}, aa*vector3{-1.0/2.0, sqrt(3.0)/2.0, 0.0}, {0.0_b, 0.0_b, lz}).spacing(aa/15.0).periodicity(2);
	
	systems::ions ions(box);
	
	ions.insert("C", {0.0_b, 0.0_b, 0.0_b});
	ions.insert("C", {0.0_b, dcc,   0.0_b});

	{
		systems::electrons electrons(env.par(), ions, box, input::config::extra_states(2), input::kpoints::grid({1, 1, 1}, false));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -11.794106663282);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         9.555702987386);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -3.535844485436);
		energy_match.check("Hartree energy",      result.energy.hartree(),       -11.139895259643);
		energy_match.check("external energy",     result.energy.external(),       15.119948124584);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.139268347279);
		energy_match.check("XC energy",           result.energy.xc(),             -4.376877703296);
		energy_match.check("XC density integral", result.energy.nvxc(),           -4.792436730841);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -19.813716465033);
		
	}
		
	fftw_cleanup();
	
	return energy_match.fail();
	
}


