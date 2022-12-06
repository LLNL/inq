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

#include <inq/inq.hpp>

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

		long nsteps = 101;
		 
		math::array<double, 1> time(nsteps);
		math::array<double, 1> dip(nsteps);
		math::array<double, 1> en(nsteps);		
	
		auto output = [&](auto data){

			auto iter = data.iter();
			
			time[iter] = data.time();
			dip[iter] = data.dipole()[0];
			en[iter] = data.energy();			

			if(iter%500 == 0 or data.last_iter()){
				auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), dip({0, iter - 1}));  
				
				std::ofstream file("spectrum.dat");
				
				for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
					file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
				}
			}
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick);

		match.check("ETRS kick: dipole step   0", dip[0],   -0.044597172494);
		match.check("ETRS kick: dipole step  10", dip[10],  -0.377148810952);
		match.check("ETRS kick: dipole step  20", dip[20],  -0.526468937415);
		match.check("ETRS kick: dipole step  30", dip[30],  -0.552455206539);
		match.check("ETRS kick: dipole step  40", dip[40],  -0.498816106782);
		match.check("ETRS kick: dipole step  50", dip[50],  -0.400067846407);
		match.check("ETRS kick: dipole step  60", dip[60],  -0.281320415502);
		match.check("ETRS kick: dipole step  70", dip[70],  -0.163199297778);
		match.check("ETRS kick: dipole step  80", dip[80],  -0.056761979779);
		match.check("ETRS kick: dipole step  90", dip[90],   0.025515738794);
		match.check("ETRS kick: dipole step 100", dip[100],  0.088155616086);		

		match.check("ETRS kick: energy step   0", en[0],   -17.563614908156);
		match.check("ETRS kick: energy step  10", en[10],  -17.563607651688);
		match.check("ETRS kick: energy step  20", en[20],  -17.563616266981);
		match.check("ETRS kick: energy step  30", en[30],  -17.563622641179);
		match.check("ETRS kick: energy step  40", en[40],  -17.563630171529);
		match.check("ETRS kick: energy step  50", en[50],  -17.563636231783);
		match.check("ETRS kick: energy step  60", en[60],  -17.563642456027);
		match.check("ETRS kick: energy step  70", en[70],  -17.563649390583);
		match.check("ETRS kick: energy step  80", en[80],  -17.563655110602);
		match.check("ETRS kick: energy step  90", en[90],  -17.563660625102);
		match.check("ETRS kick: energy step 100", en[100], -17.563665981890);	
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
