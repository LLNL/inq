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

#include <fftw3.h>

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <utils/match.hpp>
#include <operations/io.hpp>
#include <perturbations/kick.hpp>
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>

int main(int argc, char ** argv){

	using namespace inq::input;
	using inq::math::vec3d;
	
	boost::mpi3::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	inq::utils::match match(2.5e-4);

	std::vector<atom> geo;

	geo.push_back( "O" | vec3d( 0.0,      -0.553586, 0.0));
	geo.push_back( "H" | vec3d( 1.429937,  0.553586, 0.0));
  geo.push_back( "H" | vec3d(-1.429937,  0.553586, 0.0));

	inq::systems::ions ions(cell::cubic(12.0, 11.0, 10.0) | cell::finite(), geo);

  auto scf_options = scf::conjugate_gradient() | scf::energy_tolerance(1.0e-5) | scf::density_mixing() | scf::broyden_mixing();	
  
  config conf;
  
	inq::systems::electrons electrons(comm_world, ions, basis::cutoff_energy(30.0), conf);

	inq::ground_state::initialize(ions, electrons);
	auto result = inq::ground_state::calculate(ions, electrons, interaction::dft(), scf_options);
	
	match.check("total energy",        result.energy.total(),       -25.637012688764);
  match.check("kinetic energy",      result.energy.kinetic(),      12.055655438508);
  match.check("eigenvalues",         result.energy.eigenvalues,    -4.066598396189);
  match.check("Hartree energy",      result.energy.hartree,        21.255019237963);
  match.check("external energy",     result.energy.external,      -50.404984581658);
  match.check("non-local energy",    result.energy.nonlocal,       -2.732670183975);
  match.check("XC energy",           result.energy.xc,             -4.762296148152);
  match.check("XC density integral", result.energy.nvxc,           -5.494637544989);
  match.check("HF exchange energy",  result.energy.hf_exchange,     0.000000000000);
  match.check("ion-ion energy",      result.energy.ion,            -1.047736451449);

	match.check("dipole x", result.dipole[0], -0.000301716434);
	match.check("dipole y", result.dipole[1], -0.724131351204);
	match.check("dipole z", result.dipole[2], -0.000091128328);
	
	inq::operations::io::save("h2o_restart", electrons.phi_);

	fftw_cleanup(); //required for valgrid
	
	return match.fail();

}
