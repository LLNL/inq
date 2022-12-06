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

#include <fftw3.h>

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/parse_xyz.hpp>
#include <utils/match.hpp>
#include <operations/io.hpp>
#include <perturbations/kick.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>

#include <input/environment.hpp>

#include <utils/profiling.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);

	utils::match match(1e-5);

	auto box = systems::box::orthorhombic(12.0_b, 11.0_b, 10.0_b).finite().cutoff_energy(30.0_Ha);
	
	systems::ions ions(box);

	ions.insert(input::parse_xyz(config::path::unit_tests_data() + "water.xyz"));

	input::config conf;

	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	systems::electrons electrons(env.par().states(parstates), ions, box, conf);

	// Propagation without perturbation
	{
		electrons.load("h2o_restart");

		std::vector<double> energy;
		auto output = [&energy](auto data){
			energy.push_back(data.energy());
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));
		
		match.check("ETRS: energy step   0", energy[0],   -17.604152928271);
		match.check("ETRS: energy step  10", energy[10],  -17.604152928272);
		match.check("ETRS: energy step  20", energy[20],  -17.604152928272);
		match.check("ETRS: energy step  30", energy[30],  -17.604152928271);
	}

	// Propagation without perturbation
	{
		electrons.load("h2o_restart");

		std::vector<double> energy;
		auto output = [&energy](auto data){
			energy.push_back(data.energy());
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime) | input::rt::crank_nicolson());
		
		match.check("CN: energy step   0", energy[0],   -17.604152928271);
		match.check("CN: energy step  10", energy[10],  -17.604152928278);
		match.check("CN: energy step  20", energy[20],  -17.604152928294);
		match.check("CN: energy step  30", energy[30],  -17.604152928302);
	}
	
	{
		electrons.load("h2o_restart");

		auto kick = perturbations::kick{{0.1, 0.0, 0.0}};
		
		auto dipole_file = std::ofstream("dipole_etrs.dat");
		auto output = [&](auto data){
			dipole_file << data.time() << '\t' << data.dipole() << std::endl;
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick);
	}
	
	{
		electrons.load("h2o_restart");

		auto kick = perturbations::kick{{0.1, 0.0, 0.0}};

		auto dipole_file = std::ofstream("dipole_cn.dat");
		auto output = [&](auto data){
			dipole_file << data.time() << '\t' << data.dipole() << std::endl;
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime) | input::rt::crank_nicolson(), ions::propagator::fixed{}, kick);
	}
	
	fftw_cleanup(); //required for valgrind
	
	return match.fail();

}
