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

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);

	utils::match energy_match(3.0e-5);

  input::poscar vasp_file(config::path::unit_tests_data() + "bn.poscar");

	auto box = systems::box::lattice(1.0_b*vasp_file.lattice_vectors()[0], 1.0_b*vasp_file.lattice_vectors()[1], 1.0_b*vasp_file.lattice_vectors()[2]).cutoff_energy(35.0_Ha);
  
	systems::ions ions(box);
	
	ions.insert(vasp_file.geo());
  
  systems::electrons electrons(env.par(), ions, box, input::config::extra_states(3), input::kpoints::grid({2, 2, 2}, true));
	
  ground_state::initial_guess(ions, electrons);
	
  auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
	
  energy_match.check("total energy",        result.energy.total(),         -13.415825003082);
  energy_match.check("kinetic energy",      result.energy.kinetic(),         9.561515088876);
  energy_match.check("eigenvalues",         result.energy.eigenvalues,      -0.946691092608);
  energy_match.check("Hartree energy",      result.energy.hartree,           1.768136541429);
  energy_match.check("external energy",     result.energy.external,         -7.925453096022);
  energy_match.check("non-local energy",    result.energy.nonlocal,         -1.119772742664);
  energy_match.check("XC energy",           result.energy.xc,               -4.410188375663);
  energy_match.check("XC density integral", result.energy.nvxc,             -4.999253425658);
  energy_match.check("ion-ion energy",      result.energy.ion,             -11.290062419039);
	
	return energy_match.fail();
	
}


