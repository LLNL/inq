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
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>

#include <input/environment.hpp>

#include <caliper/cali.h>

int main(int argc, char ** argv){

	CALI_CXX_MARK_FUNCTION;

	using namespace inq;
	
	input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match match(1e-5);

	std::vector<input::atom> geo;

	geo.push_back( "O" | math::vec3d( 0.0,      -0.553586, 0.0));
	geo.push_back( "H" | math::vec3d( 1.429937,  0.553586, 0.0));
  geo.push_back( "H" | math::vec3d(-1.429937,  0.553586, 0.0));

	systems::ions ions(input::cell::cubic(12.0, 11.0, 10.0) | input::cell::finite(), geo);

  input::config conf;
  
  systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(30.0), conf);

	// Propagation without perturbation
	{
		operations::io::load("h2o_restart", electrons.phi_);
		
		auto result = real_time::propagate(ions, electrons, input::interaction::dft(), input::rt::num_steps(100) | input::rt::dt(0.055));
		
		match.check("energy step   0", result.energy[0],   -25.637012688816);
		match.check("energy step  10", result.energy[10],  -25.637012688717);
		match.check("energy step  20", result.energy[20],  -25.637012688605);
		match.check("energy step  30", result.energy[30],  -25.637012688485);
		match.check("energy step  40", result.energy[40],  -25.637012688359);
		match.check("energy step  50", result.energy[50],  -25.637012688229);
		match.check("energy step  60", result.energy[60],  -25.637012688095);
		match.check("energy step  70", result.energy[70],  -25.637012687958);
		match.check("energy step  80", result.energy[80],  -25.637012687819);
		match.check("energy step  90", result.energy[90],  -25.637012687679);
		match.check("energy step 100", result.energy[100], -25.637012687536);
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
		match.check("dipole x step  40", result.dipole[40][0],  -0.00035331);
		match.check("dipole y step  40", result.dipole[40][1],  -2.812590006079);
		match.check("dipole z step  40", result.dipole[40][2],  -0.00064554);
		match.check("dipole x step  50", result.dipole[50][0],  -0.000353744);
		match.check("dipole y step  50", result.dipole[50][1],  -2.812591796841);
		match.check("dipole z step  50", result.dipole[50][2],  -0.000645826);
		match.check("dipole x step  60", result.dipole[60][0],  -0.00035331);
		match.check("dipole y step  60", result.dipole[60][1],  -2.812590006079);
		match.check("dipole z step  60", result.dipole[60][2],  -0.00064554);
		match.check("dipole x step  70", result.dipole[70][0],  -0.000353744);
		match.check("dipole y step  70", result.dipole[70][1],  -2.812591042731);
		match.check("dipole z step  70", result.dipole[70][2],  -0.000645826);
		match.check("dipole x step  80", result.dipole[80][0],  -0.00035331);
		match.check("dipole y step  80", result.dipole[80][1],  -2.812591796841);
		match.check("dipole z step  80", result.dipole[80][2],  -0.00064554);
		match.check("dipole x step  90", result.dipole[90][0],  -0.00035331);
		match.check("dipole y step  90", result.dipole[90][1],  -2.812592233053);
		match.check("dipole z step  90", result.dipole[90][2],  -0.00064554);
		match.check("dipole x step 100", result.dipole[100][0], -0.000353744);
		match.check("dipole y step 100", result.dipole[100][1], -2.812592836478);
		match.check("dipole z step 100", result.dipole[100][2], -0.000645826);
		*/
	}
	
	{
		operations::io::load("h2o_restart", electrons.phi_);
		
		perturbations::kick({0.1, 0.0, 0.0}, electrons.phi_);
		
		auto result = real_time::propagate(ions, electrons, input::interaction::dft(), input::rt::num_steps(100) | input::rt::dt(0.055));
		
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

	fftw_cleanup(); //required for valgrid
	
	return match.fail();

}
