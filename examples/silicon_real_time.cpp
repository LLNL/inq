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

	auto box = systems::box::cubic(a).spacing(a/24);
	
	systems::ions ions(box);
	
	ions.insert("Si", {0.0_crys,  0.0_crys,  0.0_crys });
	ions.insert("Si", {0.25_crys, 0.25_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.5_crys,  0.0_crys });
	ions.insert("Si", {0.75_crys, 0.75_crys, 0.25_crys});
	ions.insert("Si", {0.5_crys,  0.0_crys,  0.5_crys });
	ions.insert("Si", {0.75_crys, 0.25_crys, 0.75_crys});
	ions.insert("Si", {0.0_crys,  0.5_crys,  0.5_crys });
	ions.insert("Si", {0.25_crys, 0.75_crys, 0.75_crys});

	auto nk = 2;
	
	systems::electrons electrons(env.par(), ions, box, input::kpoints::grid({nk, nk, nk}, true));

	auto functional = input::interaction::pbe();
	
	if(not electrons.try_load("silicon_restart")){
		ground_state::initial_guess(ions, electrons);
		ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::energy_tolerance(1e-4_Ha));
		ground_state::calculate(ions, electrons, functional, inq::input::scf::energy_tolerance(1e-8_Ha));
		electrons.save("silicon_restart");
	}

	auto kick = perturbations::kick{box.cell(), {0.01, 0.0, 0.0}, perturbations::gauge::velocity};
	
	auto const dt = 0.065;
	long nsteps = 413.41373/dt;
	
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
	
	real_time::propagate<>(ions, electrons, output, functional, input::rt::num_steps(nsteps) | input::rt::dt(dt*1.0_atomictime) | input::rt::etrs(), ions::propagator::fixed{}, kick);
	
	return energy_match.fail();
	
}


