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
		
		auto result = real_time::propagate<>(ions, electrons, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));
		
		match.check("ETRS: energy step   0", result.energy[0],   -25.637012688816);
		match.check("ETRS: energy step  10", result.energy[10],  -25.637012688717);
		match.check("ETRS: energy step  20", result.energy[20],  -25.637012688605);
		match.check("ETRS: energy step  30", result.energy[30],  -25.637012688485);
	}

	// Propagation without perturbation
	{
		electrons.load("h2o_restart");
		
		auto result = real_time::propagate<>(ions, electrons, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime) | input::rt::crank_nicolson());
		
		match.check("CN: energy step   0", result.energy[0],   -25.637012688816);
		match.check("CN: energy step  10", result.energy[10],  -25.637012688717);
		match.check("CN: energy step  20", result.energy[20],  -25.637012688605);
		match.check("CN: energy step  30", result.energy[30],  -25.637012688485);
	}
	
	{
		electrons.load("h2o_restart");

		for(auto phi : electrons.lot()) perturbations::kick({0.1, 0.0, 0.0}, phi.fields());
		
		auto result = real_time::propagate<>(ions, electrons, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));
		
		{
			auto dipole_file = std::ofstream("dipole_etrs.dat");
			
			for(unsigned ii = 0; ii < result.dipole.size(); ii++){
				dipole_file << result.time[ii] << '\t' << result.dipole[ii] << std::endl;
			}
		}
	}
	
	{
		electrons.load("h2o_restart");
		
		for(auto phi : electrons.lot()) perturbations::kick({0.1, 0.0, 0.0}, phi.fields());
		
		auto result = real_time::propagate<>(ions, electrons, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime) | input::rt::crank_nicolson());
		
		{
			auto dipole_file = std::ofstream("dipole_cn.dat");
			
			for(unsigned ii = 0; ii < result.dipole.size(); ii++){
				dipole_file << result.time[ii] << '\t' << result.dipole[ii] << std::endl;
			}
		}
	}
	
	fftw_cleanup(); //required for valgrind
	
	return match.fail();

}
