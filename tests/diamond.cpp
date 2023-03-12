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

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <observables/kinetic_energy_density.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;

	input::environment env(argc, argv);

	utils::match energy_match(3.0e-5);

	std::vector<input::atom> geo;

	auto a =  3.567095_A;

	auto box = systems::box::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}).cutoff_energy(35.0_Ha);
	
	systems::ions ions(box);
	
	ions.insert("C", {0.0_crys,  0.0_crys,  0.0_crys });
	ions.insert("C", {0.25_crys, 0.25_crys, 0.25_crys});

	{
		systems::electrons electrons(env.par(), ions, box, input::config::extra_states(3), input::kpoints::grid({1, 1, 1}, false));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -10.949196617732);
		energy_match.check("kinetic energy",      result.energy.kinetic(),        11.411454690719);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),     0.795474829167);
		energy_match.check("Hartree energy",      result.energy.hartree(),         1.473883973621);
		energy_match.check("external energy",     result.energy.external(),       -7.034444637769);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.496762367780);
		energy_match.check("XC energy",           result.energy.xc(),             -4.568604147920);
		energy_match.check("XC density integral", result.energy.nvxc(),           -5.032540803244);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
		
		auto ked = observables::kinetic_energy_density(electrons);
		
		energy_match.check("kinetic energy", operations::integral(ked), 11.411455188639);
	}

	systems::electrons electrons(env.par(), ions, box, input::config::extra_states(3), input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);

	{
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));

		energy_match.check("total energy",        result.energy.total(),         -12.041228146904);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.513350495571);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.761571208702);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.974854048324);
		energy_match.check("external energy",     result.energy.external(),       -5.856397438020);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.593446844490);
		energy_match.check("XC energy",           result.energy.xc(),             -4.344864279688);
		energy_match.check("XC density integral", result.energy.nvxc(),           -4.774785518412);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),     0.000000000000);		
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
		
		auto ked = observables::kinetic_energy_density(electrons);
		
		energy_match.check("kinetic energy", operations::integral(ked), 8.513350448818);
	}

	{
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe0(), inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(0));
		
		energy_match.check("total energy",        result.energy.total(),         -11.569230679378);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.513350460378);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.715154057693);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.974854033969);
		energy_match.check("external energy",     result.energy.external(),       -5.856397409990);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.593446829105);
		energy_match.check("XC energy",           result.energy.xc(),             -3.357693539830);
		energy_match.check("XC density integral", result.energy.nvxc(),           -3.698021814519);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),    -0.515173266198);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}
	
	{
		auto result = ground_state::calculate(ions, electrons, input::interaction::pbe0(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -11.571610770029);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.386909236299);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.752041292771);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.939332482432);
		energy_match.check("external energy",     result.energy.external(),       -5.736703986825);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.566028479461);
		energy_match.check("XC energy",           result.energy.xc(),             -3.345210887801);
		energy_match.check("XC density integral", result.energy.nvxc(),           -3.684513015509);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),    -0.515185006069);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}

	{
		auto result = ground_state::calculate(ions, electrons, input::interaction::hartree_fock(), inq::input::scf::steepest_descent() | inq::input::scf::energy_tolerance(1e-8_Ha));

		energy_match.check("total energy",        result.energy.total(),          -9.788709725748);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.151819376871);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.249826944617);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.880126063939);
		energy_match.check("external energy",     result.energy.external(),       -5.499063363813);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.510900262733);
		energy_match.check("XC energy",           result.energy.xc(),              0.000000000000);
		energy_match.check("XC density integral", result.energy.nvxc(),            0.000000000000);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),    -2.075967411411);		
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}
	
	fftw_cleanup();
	
	return energy_match.fail();
	
}


