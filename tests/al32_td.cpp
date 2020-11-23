/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa

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
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <input/parse_xyz.hpp>
#include <config/path.hpp>

#include <real_time/propagate.hpp>

#include<fstream>

int main(int argc, char ** argv){

	CALI_CXX_MARK_FUNCTION;

	using namespace inq;
	
	input::environment env(argc, argv);

	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(4.0e-6);

	std::vector<input::atom> geo = input::parse_xyz(config::path::unit_tests_data() + "al32.xyz");

	geo.emplace_back("H" | math::vec3d(0.00000, 1.91325, 1.91325));

	systems::ions ions(input::cell::cubic(2*7.6524459), geo);
	
	input::config conf;
	
	conf.excess_charge = -1;
	conf.extra_states = 8;
	conf.temperature = 0.00095004347; //300 K
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0), conf);
	
	ground_state::initialize(ions, electrons);

	auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(200));
	
	energy_match.check("total energy",        result.energy.total(),        -73.533747057250);
	energy_match.check("kinetic energy",      result.energy.kinetic(),       28.023995211588);
	energy_match.check("eigenvalues",         result.energy.eigenvalues,      6.501664582729);
	energy_match.check("Hartree energy",      result.energy.hartree,          0.122273970540);
	energy_match.check("external energy",     result.energy.external,         1.593742691403);
	energy_match.check("non-local energy",    result.energy.nonlocal,        11.261163563655);
	energy_match.check("XC energy",           result.energy.xc,             -34.568675453999);
	energy_match.check("XC density integral", result.energy.nvxc,           -34.621784824997);
	energy_match.check("ion-ion energy",      result.energy.ion,            -79.966247040437);

	inq::operations::io::save("al32_restart", electrons.phi_);

	auto dt = 0.055;

	ions.velocities()[ions.geo().num_atoms() - 1] = math::vec3d(0.1, 0.0, 0.0);

	{
		auto propagation = real_time::propagate(
			ions, electrons, 
			input::interaction::non_interacting(), input::rt::num_steps(1000) | input::rt::dt(dt), 
			real_time::impulsive_ions{}
		);

		auto ofs = std::ofstream{"al32_v0.1.dat"}; ofs<<"# distance (au), energy (au)\n";

		for(std::size_t i = 0; i != propagation.ions.size(); ++i)
			ofs<< propagation.ions[i].geo().coordinates()[ions.geo().num_atoms()-1][0] <<'\t'<< propagation.energy[i] <<'\n';
	}
	fftw_cleanup(); //required for valgrid
	
	return energy_match.fail();
	
}

