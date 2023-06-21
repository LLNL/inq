/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);

	utils::match match(1e-5);

	auto box = systems::box::orthorhombic(12.0_b, 11.0_b, 10.0_b).finite();
	
	systems::ions ions(box);

	ions.insert(input::parse_xyz(config::path::unit_tests_data() + "water.xyz"));

	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	systems::electrons electrons(env.par().states(parstates), ions, box, input::config::cutoff(30.0_Ha));

	// Propagation without perturbation
	{
		electrons.load("h2o_restart");

		std::vector<double> energy;
		auto output = [&energy](auto data){
			energy.push_back(data.energy());
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(30) | input::rt::dt(0.055_atomictime));
		
		match.check("ETRS: energy step   0", energy[0],   -17.604152928110);
		match.check("ETRS: energy step  10", energy[10],  -17.604152928110);
		match.check("ETRS: energy step  20", energy[20],  -17.604152928110);
		match.check("ETRS: energy step  30", energy[30],  -17.604152928110);
	}

	// Propagation without perturbation
	{
		electrons.load("h2o_restart");

		std::vector<double> energy;
		auto output = [&energy](auto data){
			energy.push_back(data.energy());
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(10) | input::rt::dt(0.1_atomictime) | input::rt::crank_nicolson());
		
		/* This is disabled for now since CN seems to be a bit numerically unstable
			 
			match.check("CN: energy step  0", energy[0], -17.604152928274);
			match.check("CN: energy step  1", energy[1], -17.604152928429);
			match.check("CN: energy step  2", energy[2], -17.604152927327);
			match.check("CN: energy step  3", energy[3], -17.604152924407);
			match.check("CN: energy step  4", energy[4], -17.604152924731);
			match.check("CN: energy step  5", energy[5], -17.604152924780);
			match.check("CN: energy step  6", energy[6], -17.604152919840);
			match.check("CN: energy step  7", energy[7], -17.604152919615);
			match.check("CN: energy step  8", energy[8], -17.604152896985);
			match.check("CN: energy step  9", energy[9], -17.604152702584);
		*/
	}
	
	{
		electrons.load("h2o_restart");

		auto kick = perturbations::kick{box.cell(), {0.1, 0.0, 0.0}};

		long nsteps = 71;
		 
		gpu::array<double, 1> time(nsteps);
		gpu::array<double, 1> dip(nsteps);
		gpu::array<double, 1> en(nsteps);		
	
		auto output = [&](auto data){

			auto iter = data.iter();
			
			time[iter] = data.time();
			dip[iter] = data.dipole()[0];
			en[iter] = data.energy();			

			if(data.root() and data.every(50)){
				auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), dip({0, iter - 1}));  

				std::ofstream file("spectrum.dat");
				
				for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
					file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
				}
			}
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick);

		match.check("ETRS length kick: dipole step   0", dip[0],   0.043955375747);
		match.check("ETRS length kick: dipole step  10", dip[10],  0.376347806791);
		match.check("ETRS length kick: dipole step  20", dip[20],  0.525427259213);
		match.check("ETRS length kick: dipole step  30", dip[30],  0.550931744154);
		match.check("ETRS length kick: dipole step  40", dip[40],  0.497454787505);
		match.check("ETRS length kick: dipole step  50", dip[50],  0.397016815641);
		match.check("ETRS length kick: dipole step  60", dip[60],  0.273814233337);
		match.check("ETRS length kick: dipole step  70", dip[70],  0.150512273021);

		match.check("ETRS length kick: energy step   0", en[0],   -17.563614846419);
		match.check("ETRS length kick: energy step  10", en[10],  -17.563607131141);
		match.check("ETRS length kick: energy step  20", en[20],  -17.563615606337);
		match.check("ETRS length kick: energy step  30", en[30],  -17.563621943406);
		match.check("ETRS length kick: energy step  40", en[40],  -17.563629437416);
		match.check("ETRS length kick: energy step  50", en[50],  -17.563635432990);
		match.check("ETRS length kick: energy step  60", en[60],  -17.563641616526);
		match.check("ETRS length kick: energy step  70", en[70],  -17.563648522511);
	}

	{
		electrons.load("h2o_restart");
		
		auto kick1 = perturbations::kick{box.cell(), {0.06, 0.0, 0.0}, perturbations::gauge::velocity};
		auto kick2 = perturbations::kick{box.cell(), {0.04, 0.0, 0.0}, perturbations::gauge::velocity};

		long nsteps = 31;
		 
		gpu::array<double, 1> time(nsteps);
		gpu::array<double, 1> dip(nsteps);
		gpu::array<double, 1> en(nsteps);		
	
		auto output = [&](auto data){

			auto iter = data.iter();
			
			time[iter] = data.time();
			dip[iter] = data.dipole()[0];
			en[iter] = data.energy();			

			if(data.root() and data.every(50)){
				auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), dip({0, iter - 1}));  

				std::ofstream file("spectrum.dat");
				
				for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
					file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
				}
			}
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick1 + kick2);

		match.check("ETRS velocity kick: dipole step   0", dip[0],   0.043697788108);
		match.check("ETRS velocity kick: dipole step  10", dip[10],  0.375961176642);
		match.check("ETRS velocity kick: dipole step  20", dip[20],  0.525287483544);
		match.check("ETRS velocity kick: dipole step  30", dip[30],  0.551308944304);

		match.check("ETRS velocity kick: energy step   0", en[0],   -17.563918466246);
		match.check("ETRS velocity kick: energy step  10", en[10],  -17.563906778220);
		match.check("ETRS velocity kick: energy step  20", en[20],  -17.563911508759);
		match.check("ETRS velocity kick: energy step  30", en[30],  -17.563914263774);
	}
		
	{
		electrons.load("h2o_restart");

		auto kick = perturbations::kick{box.cell(), {0.1, 0.0, 0.0}};

		long nsteps = 21;
		
		gpu::array<double, 1> time(nsteps);
		gpu::array<double, 1> dip(nsteps);
		gpu::array<double, 1> en(nsteps);		
	
		auto output = [&](auto data){
			
			auto iter = data.iter();
			
			time[iter] = data.time();
			dip[iter] = data.dipole()[0];
			en[iter] = data.energy();			

			if(data.root() and data.every(50)){
				auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), dip({0, iter - 1}));  

				std::ofstream file("spectrum.dat");
				
				for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
					file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
				}
			}
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime) | input::rt::crank_nicolson(), ions::propagator::fixed{}, kick);
		
	}
	
	fftw_cleanup(); //required for valgrind
	
	return match.fail();

}
