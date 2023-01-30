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

	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	systems::electrons electrons(env.par().states(parstates), ions, box);

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
		
		match.check("CN: energy step  0", energy[0], -17.604152928110);
		match.check("CN: energy step  1", energy[1], -17.604152928110);
		match.check("CN: energy step  2", energy[2], -17.604152928109);
		match.check("CN: energy step  3", energy[3], -17.604152928089);
		match.check("CN: energy step  4", energy[4], -17.604152928086);
		match.check("CN: energy step  5", energy[5], -17.604152928082);
		match.check("CN: energy step  6", energy[6], -17.604152928081);
		match.check("CN: energy step  7", energy[7], -17.604152928076);
		match.check("CN: energy step  8", energy[8], -17.604152928064);
		match.check("CN: energy step  9", energy[9], -17.604152928054);		
	}
	
	{
		electrons.load("h2o_restart");

		auto kick = perturbations::kick{box.cell(), {0.1, 0.0, 0.0}};

		long nsteps = 71;
		 
		math::array<double, 1> time(nsteps);
		math::array<double, 1> dip(nsteps);
		math::array<double, 1> en(nsteps);		
	
		auto output = [&](auto data){

			auto iter = data.iter();
			
			time[iter] = data.time();
			dip[iter] = data.dipole()[0];
			en[iter] = data.energy();			

			if(data.every(50)){
				auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), dip({0, iter - 1}));  

				std::ofstream file("spectrum.dat");
				
				for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
					file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
				}
			}
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick);

		match.check("ETRS length kick: dipole step   0", dip[0],    0.043955102932);
		match.check("ETRS length kick: dipole step  10", dip[10],   0.376347925882);
		match.check("ETRS length kick: dipole step  20", dip[20],   0.525428390792);
		match.check("ETRS length kick: dipole step  30", dip[30],   0.550931444557);
		match.check("ETRS length kick: dipole step  40", dip[40],   0.497454742178);
		match.check("ETRS length kick: dipole step  50", dip[50],   0.397014962050);
		match.check("ETRS length kick: dipole step  60", dip[60],   0.273822734950);
		match.check("ETRS length kick: dipole step  70", dip[70],   0.150516269346);

		match.check("ETRS length kick: energy step   0", en[0],   -17.563615649009);
		match.check("ETRS length kick: energy step  10", en[10],  -17.563615382735);
		match.check("ETRS length kick: energy step  20", en[20],  -17.563630261556);
		match.check("ETRS length kick: energy step  30", en[30],  -17.563640944592);
		match.check("ETRS length kick: energy step  40", en[40],  -17.563653420614);
		match.check("ETRS length kick: energy step  50", en[50],  -17.563663830708);
		match.check("ETRS length kick: energy step  60", en[60],  -17.563673934214);
		match.check("ETRS length kick: energy step  70", en[70],  -17.563684356007);
	}

	{
		electrons.load("h2o_restart");
		
		auto kick1 = perturbations::kick{box.cell(), {0.06, 0.0, 0.0}, perturbations::gauge::velocity};
		auto kick2 = perturbations::kick{box.cell(), {0.04, 0.0, 0.0}, perturbations::gauge::velocity};

		long nsteps = 31;
		 
		math::array<double, 1> time(nsteps);
		math::array<double, 1> dip(nsteps);
		math::array<double, 1> en(nsteps);		
	
		auto output = [&](auto data){

			auto iter = data.iter();
			
			time[iter] = data.time();
			dip[iter] = data.dipole()[0];
			en[iter] = data.energy();			

			if(data.every(50)){
				auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), dip({0, iter - 1}));  

				std::ofstream file("spectrum.dat");
				
				for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
					file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
				}
			}
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick1 + kick2);

		match.check("ETRS velocity kick: dipole step   0", dip[0],   0.043697635032);
		match.check("ETRS velocity kick: dipole step  10", dip[10],  0.375961121545);
		match.check("ETRS velocity kick: dipole step  20", dip[20],  0.525287647225);
		match.check("ETRS velocity kick: dipole step  30", dip[30],  0.551308998120);

		match.check("ETRS velocity kick: energy step   0", en[0],   -17.563918857038);
		match.check("ETRS velocity kick: energy step  10", en[10],  -17.563910876489);
		match.check("ETRS velocity kick: energy step  20", en[20],  -17.563916148234);
		match.check("ETRS velocity kick: energy step  30", en[30],  -17.563918743496);
	}
		
	{
		electrons.load("h2o_restart");

		auto kick = perturbations::kick{box.cell(), {0.1, 0.0, 0.0}};

		auto dipole_file = std::ofstream("dipole_cn.dat");
		auto output = [&](auto data){
			dipole_file << data.time() << '\t' << data.dipole() << std::endl;
		};
		
		real_time::propagate<>(ions, electrons, output, input::interaction::lda(), input::rt::num_steps(10) | input::rt::dt(0.1_atomictime) | input::rt::crank_nicolson(), ions::propagator::fixed{}, kick);
	}
	
	fftw_cleanup(); //required for valgrind
	
	return match.fail();

}
