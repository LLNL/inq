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
#include <input/atom.hpp>
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

	std::vector<input::atom> geo;

	geo.push_back( "O" | math::vector3<double>( 0.0,      -0.553586, 0.0));
	geo.push_back( "H" | math::vector3<double>( 1.429937,  0.553586, 0.0));
	geo.push_back( "H" | math::vector3<double>(-1.429937,  0.553586, 0.0));

	auto box = systems::box::orthorhombic(12.0_b, 11.0_b, 10.0_b).finite().cutoff_energy(30.0_Ha);
	
	systems::ions ions(box, geo);

	input::config conf;

	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	systems::electrons electrons(env.par().states(parstates), ions, box, conf);

	// Propagation without perturbation
	{
		electrons.load("h2o_restart");
		
		auto result = real_time::propagate<>(ions, electrons, input::interaction::dft(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));
		
		match.check("energy step   0", result.energy[0],   -25.637012688816);
		match.check("energy step  10", result.energy[10],  -25.637012688717);
		match.check("energy step  20", result.energy[20],  -25.637012688605);
		match.check("energy step  30", result.energy[30],  -25.637012688485);

		/*
		match.check("dipole x step   0", result.dipole[0][0],   -0.00035331);
		match.check("dipole y step   0", result.dipole[0][1],   -2.812602083171);
		match.check("dipole z step   0", result.dipole[0][2],   -0.00064554);
		match.check("dipole x step  10", result.dipole[10][0],  -0.000353744);
		match.check("dipole y step  10", result.dipole[10][1],  -2.812598839225);
		match.check("dipole z step  10", result.dipole[10][2],  -0.000645826);
		match.check("dipole x step  20", result.dipole[20][0],  -0.00035331);
		match.check("dipole y step  20", result.dipole[20][1],  -2.812591107786);
		match.check("dipole z step  20", result.dipole[20][2],  -0.00064554);
		match.check("dipole x step  30", result.dipole[30][0],  -0.000353744);
		match.check("dipole y step  30", result.dipole[30][1],  -2.812589627003);
		match.check("dipole z step  30", result.dipole[30][2],  -0.000645826);
		*/
	}
	
	{
		electrons.load("h2o_restart");

		for(auto phi : electrons.lot()) perturbations::kick({0.1, 0.0, 0.0}, phi.fields());
		
		auto result = real_time::propagate<>(ions, electrons, input::interaction::dft(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));
		
		/*		match.check("energy step  0", result.energy[0], -16.903925978590);
		match.check("energy step 10", result.energy[10], -16.904635586794);
		*/
		{
			auto dipole_file = std::ofstream("dipole.dat");
			
			for(unsigned ii = 0; ii < result.dipole.size(); ii++){
				dipole_file << result.time[ii] << '\t' << result.dipole[ii] << std::endl;
			}
		}
	}

	fftw_cleanup(); //required for valgrind
	
	return match.fail();

}
