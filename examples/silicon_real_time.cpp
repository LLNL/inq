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

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);

	utils::match energy_match(3.0e-5);

	std::vector<input::atom> geo;

	auto a = 10.18_b;

	auto box = systems::box::cubic(a).cutoff_energy(30.0_Ha);
	
	systems::ions ions(box);
	
	ions.insert("Si", {0.0_crys,  0.0_crys,  0.0_crys });
	ions.insert("Si", {0.25_crys, 0.25_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.5_crys,  0.0_crys });
	ions.insert("Si", {0.75_crys, 0.75_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.0_crys,  0.5_crys });
	ions.insert("Si", {0.75_crys, 0.25_crys, 0.75_crys});
	ions.insert("Si", {0.0_crys,  0.5_crys,  0.5_crys });
	ions.insert("Si", {0.25_crys, 0.75_crys, 0.75_crys});

	int kpoint_par = 1;
	if(env.par().size()%2 == 0) kpoint_par = 2;
	
	systems::electrons electrons(env.par().kpoints(kpoint_par), ions, box, input::kpoints::grid({4, 4, 4}, true));

	if(not electrons.try_load("silicon_restart")){
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::energy_tolerance(1e-8_Ha));
		electrons.save("silicon_restart");
	}

	auto kick = perturbations::kick{box.cell(), {0.01, 0.0, 0.0}, perturbations::gauge::velocity};
	
	long nsteps = 10001;
	
	math::array<double, 1> time(nsteps);
	math::array<double, 1> cur(nsteps);
	math::array<double, 1> en(nsteps);		

	std::ofstream file;
	if(electrons.root()) file.open("current.dat");
	
	auto output = [&](auto data){
		
		auto iter = data.iter();
		
		time[iter] = data.time();
		cur[iter] = data.current()[0];
		en[iter] = data.energy();

		if(data.root()) file << time[iter] << '\t' << cur[iter] << std::endl;
		
		if(data.root() and data.every(50)){
			auto spectrum = observables::spectrum(20.0_eV, 0.01_eV, time({0, iter - 1}), cur({0, iter - 1}));  
			
			std::ofstream file("spectrum.dat");
			
			for(int ifreq = 0; ifreq < spectrum.size(); ifreq++){
				file << ifreq*in_atomic_units(0.01_eV) << '\t' << real(spectrum[ifreq]) << '\t' << imag(spectrum[ifreq]) << std::endl;
			}
		}
	};
	
	real_time::propagate<>(ions, electrons, output, input::interaction::pbe(), input::rt::num_steps(nsteps) | input::rt::dt(0.055_atomictime), ions::propagator::fixed{}, kick);
	
	return energy_match.fail();
	
}


