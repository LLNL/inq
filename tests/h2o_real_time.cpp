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
#include <utils/match.hpp>
#include <operations/io.hpp>
#include <perturbations/kick.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>

int main(int argc, char ** argv){

	boost::mpi3::environment env(argc, argv);

	utils::match energy_match(2.5e-4);

	std::vector<input::atom> geo;

	geo.push_back( "O" | math::vec3d( 0.0,      -0.553586, 0.0));
	geo.push_back( "H" | math::vec3d( 1.429937,  0.553586, 0.0));
  geo.push_back( "H" | math::vec3d(-1.429937,  0.553586, 0.0));

	systems::ions ions(input::cell::cubic(12.0, 11.0, 10.0) | input::cell::finite(), geo);

  auto scf_options = input::scf::conjugate_gradient() | input::scf::energy_tolerance(1.0e-5) | input::scf::density_mixing();
  
  input::config conf;
  
  systems::electrons electrons(ions, input::basis::cutoff_energy(20.0), conf);

  auto energy = ground_state::calculate(electrons, input::interaction::dft(), scf_options);
  
  energy_match.check("total energy",        energy.total(),         -25.885010460377);
  energy_match.check("kinetic energy",      energy.kinetic(),        12.011756577795);
  energy_match.check("eigenvalues",         energy.eigenvalues,      -4.247001211716);
  energy_match.check("Hartree energy",      energy.hartree,          21.211005609920);
  energy_match.check("external energy",     energy.external,        -50.383959864425);
  energy_match.check("non-local energy",    energy.nonlocal,         -2.809830867482);
  energy_match.check("XC energy",           energy.xc,               -4.866245464735);
  energy_match.check("XC density integral", energy.nvxc,             -5.486978277443);
  energy_match.check("HF exchange energy",  energy.hf_exchange,       0.0);
  energy_match.check("ion-ion energy",      energy.ion,              -1.047736451449);

	operations::io::save("h2o_restart", electrons.phi_);

	operations::io::load("h2o_restart", electrons.phi_);

  perturbations::kick({1.0, 0.0, 0.0}, electrons.phi_);
											 
  auto result = real_time::propagate(electrons);

  energy_match.check("energy step  0", result.energy[0], -16.903925978590);
  energy_match.check("energy step 10", result.energy[10], -16.904635586794);

  {
    auto dipole_file = std::ofstream("dipole.dat");

    for(unsigned ii = 0; ii < result.dipole.size(); ii++){
      dipole_file << ii << '\t' << result.time[ii] << '\t' << result.dipole[ii] << std::endl;
    }
  }
  
	return energy_match.fail();
}
