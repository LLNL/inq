/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2023 Xavier Andrade

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
	
	utils::match match(3.0e-5);

	auto box = systems::box::cubic(15.0_b).finite().cutoff_energy(30.0_Ha);

	systems::ions ions(box);

	auto distance = 121.0_pm;
	
	ions.insert("O", {-distance/2, 0.0_b, 0.0_b});
	ions.insert("O", {distance/2, 0.0_b, 0.0_b});	

	systems::electrons electrons(env.par(), ions, box, input::config::spin_polarized() | input::config::temperature(1000.0_K) | input::config::extra_states(3));
	ground_state::initial_guess(ions, electrons);
		
	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe());
	
	match.check("total energy",        result.energy.total()    ,  -34.829848079842);
	match.check("kinetic energy",      result.energy.kinetic()  ,   21.076601720040);
	match.check("eigenvalues",         result.energy.eigenvalues,   -7.304540681975);
	match.check("external energy",     result.energy.external   ,  -99.306470054754);
	match.check("non-local energy",    result.energy.nonlocal   ,   -4.538047573763);
	match.check("XC energy",           result.energy.xc,           -10.061335310136);
	match.check("XC density integral", result.energy.nvxc,          -9.047200319683);
	match.check("HF exchange energy",  result.energy.hf_exchange,    0.000000000000);
	match.check("ion-ion energy",      result.energy.ion,           15.744115365679);

	return match.fail();
	
}

