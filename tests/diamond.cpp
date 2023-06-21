/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
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

	auto a = 3.567095_A;
	systems::ions ions(systems::cell::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}));
	
	ions.insert_fractional("C", {0.0,  0.0,  0.0 });
	ions.insert_fractional("C", {0.25, 0.25, 0.25});

	{
		systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(35.0_Ha).extra_states(3), input::kpoints::gamma() + input::kpoints::point({-0.5, -0.5, -0.5}, 0.0));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -10.897044233602);
		energy_match.check("kinetic energy",      result.energy.kinetic(),        11.213776066809);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),     0.828607402452);
		energy_match.check("Hartree energy",      result.energy.hartree(),         1.454295391985);
		energy_match.check("external energy",     result.energy.external(),       -6.932275373959);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -1.343071950581);
		energy_match.check("XC energy",           result.energy.xc(),             -4.555044239254);
		energy_match.check("XC density integral", result.energy.nvxc(),           -5.018412123788);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
		
		auto all_eigenvalues = parallel::gather(+electrons.eigenvalues().flatted(), electrons.kpin_states_part(), electrons.kpin_states_comm(), 0);

		if(electrons.kpin_states_comm().root()){
			energy_match.check("gamma homo",          all_eigenvalues[3],      0.309230017411);
			energy_match.check("gamma lumo",          all_eigenvalues[4],      0.503210619796);
			energy_match.check("0.5,0.5,0.5 homo",    all_eigenvalues[7 + 3],  0.201622641335);
			energy_match.check("0.5,0.5,0.5 lumo",    all_eigenvalues[7 + 4],  0.604793158669);
		}

		auto ked = observables::kinetic_energy_density(electrons);
		
		energy_match.check("kinetic energy", operations::integral(ked), 11.213776066809);
	}
	
	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(35.0_Ha).extra_states(3), input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);

	{
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));

		energy_match.check("total energy",        result.energy.total(),         -11.999181534999);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.334311155110);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.731834947260);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.956736892010);
		energy_match.check("external energy",     result.energy.external(),       -5.757998958933);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.464380542859);
		energy_match.check("XC energy",           result.energy.xc(),             -4.333125951724);
		energy_match.check("XC density integral", result.energy.nvxc(),           -4.757240384597);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),     0.000000000000);		
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
		
		auto ked = observables::kinetic_energy_density(electrons);
		
		energy_match.check("kinetic energy", operations::integral(ked), 8.334311155110);
	}

	if(energy_match.fail()) return energy_match.fail();		

	{
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe0(), inq::options::ground_state{}.steepest_descent().scf_steps(0));
		
		energy_match.check("total energy",        result.energy.total(),         -11.529999323732);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.334311155110);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.690589480039);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.956736892010);
		energy_match.check("external energy",     result.energy.external(),       -5.757998958933);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.464380542859);
		energy_match.check("XC energy",           result.energy.xc(),             -3.349059614094);
		energy_match.check("XC density integral", result.energy.nvxc(),           -3.686226664651);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),    -0.514884126363);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}

	if(energy_match.fail()) return energy_match.fail();

	{
		auto result = ground_state::calculate(ions, electrons, options::theory{}.hartree_fock(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));

		energy_match.check("total energy",        result.energy.total(),          -9.755173389899);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.028934085593);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.227947505918);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.867423458775);
		energy_match.check("external energy",     result.energy.external(),       -5.425380321045);
		energy_match.check("non-local energy",    result.energy.nonlocal(),       -0.416504781221);
		energy_match.check("XC energy",           result.energy.xc(),              0.000000000000);
		energy_match.check("XC density integral", result.energy.nvxc(),            0.000000000000);
		energy_match.check("HF exchange energy",  result.energy.hf_exchange(),    -2.074921703397);		
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}
	
	fftw_cleanup();
	
	return energy_match.fail();
	
}


