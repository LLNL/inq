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

	double a = 10.18;
	
	geo.push_back( "Si" | a*math::vector3<double>(0.0,  0.0,  0.0 ));
	geo.push_back( "Si" | a*math::vector3<double>(0.25, 0.25, 0.25));
	geo.push_back( "Si" | a*math::vector3<double>(0.5,  0.5,  0.0 ));
	geo.push_back( "Si" | a*math::vector3<double>(0.75, 0.75, 0.25));
	geo.push_back( "Si" | a*math::vector3<double>(0.5,  0.0,  0.5 ));
	geo.push_back( "Si" | a*math::vector3<double>(0.75, 0.25, 0.75));
	geo.push_back( "Si" | a*math::vector3<double>(0.0,  0.5,  0.5 ));
	geo.push_back( "Si" | a*math::vector3<double>(0.25, 0.75, 0.75));

	auto box = systems::box::cubic(a*1.0_b).cutoff_energy(25.0_Ha);
	
	systems::ions ions(box, geo);
	
	input::config conf;
	
	conf.extra_states = 0;
	
	systems::electrons electrons(env.dist(), ions, box, conf, input::kpoints::grid({2, 1, 1}, true));
	
	ground_state::initial_guess(ions, electrons);
	auto result = ground_state::calculate(ions, electrons, input::interaction::non_interacting(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
	
	energy_match.check("total energy",     result.energy.total()    , -23.834202307265);
	energy_match.check("kinetic energy",   result.energy.kinetic()  ,  14.428064504524);
	energy_match.check("eigenvalues",      result.energy.eigenvalues,   7.649418187835);
	energy_match.check("external energy",  result.energy.external   , -12.019317998765);
	energy_match.check("non-local energy", result.energy.nonlocal   ,   5.240671682076);
	energy_match.check("ion-ion energy",   result.energy.ion        , -31.483620495100);
	
	electrons.save("silicon_restart");

	auto ked = observables::kinetic_energy_density(electrons);

	energy_match.check("kinetic energy", operations::integral(ked), 14.428064504524);
	
	fftw_cleanup(); //required for valgrid
	
	return energy_match.fail();
	
}

