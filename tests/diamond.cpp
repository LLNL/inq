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
	
	utils::match energy_match(3.0e-5);

	auto a = 3.567095_A;
	systems::ions ions(systems::cell::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}));
	
	ions.insert_fractional("C", {0.0,  0.0,  0.0 });
	ions.insert_fractional("C", {0.25, 0.25, 0.25});

	{
		systems::electrons electrons(ions, options::electrons{}.cutoff(35.0_Ha).extra_states(3), input::kpoints::gamma() + input::kpoints::point({-0.5, -0.5, -0.5}, 0.0));
		
		ground_state::initial_guess(ions, electrons);
		
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
		
		energy_match.check("total energy",        result.energy.total(),         -10.949196617732);
		energy_match.check("kinetic energy",      result.energy.kinetic(),        11.411454690719);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),     0.795474829167);
		energy_match.check("Hartree energy",      result.energy.hartree(),         1.473883973621);
		energy_match.check("external energy",     result.energy.external(),       -7.034444637769);
		energy_match.check("non-local energy",    result.energy.non_local(),       -1.496762367780);
		energy_match.check("XC energy",           result.energy.xc(),             -4.568604147920);
		energy_match.check("XC density integral", result.energy.nvxc(),           -5.032540803244);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);

		auto all_eigenvalues = parallel::gather(+electrons.eigenvalues().flatted(), electrons.kpin_states_part(), electrons.kpin_states_comm(), 0);

		if(electrons.kpin_states_comm().root()){
			energy_match.check("gamma homo",          all_eigenvalues[3],  0.303880260047);
			energy_match.check("gamma lumo",          all_eigenvalues[4],  0.496554623841);
			energy_match.check("0.5,0.5,0.5 homo",    all_eigenvalues[7 + 3],  0.195774173500);
			energy_match.check("0.5,0.5,0.5 lumo",    all_eigenvalues[7 + 4],  0.597476824600);
		}

		auto ked = observables::kinetic_energy_density(electrons);
		
		energy_match.check("kinetic energy", operations::integral(ked), 11.411455188639);
	}
	
	systems::electrons electrons(ions, options::electrons{}.cutoff(35.0_Ha).extra_states(3), input::kpoints::grid({2, 2, 2}, true));
	
	ground_state::initial_guess(ions, electrons);

	{
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));

		energy_match.check("total energy",        result.energy.total(),         -12.041228146904);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.513350495571);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.761571208702);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.974854048324);
		energy_match.check("external energy",     result.energy.external(),       -5.856397438020);
		energy_match.check("non-local energy",    result.energy.non_local(),       -0.593446844490);
		energy_match.check("XC energy",           result.energy.xc(),             -4.344864279688);
		energy_match.check("XC density integral", result.energy.nvxc(),           -4.774785518412);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),     0.000000000000);		
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
		
		auto ked = observables::kinetic_energy_density(electrons);
		
		energy_match.check("kinetic energy", operations::integral(ked), 8.513350448818);
	}

	if(energy_match.fail()) return energy_match.fail();		

	{
		auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe0(), inq::options::ground_state{}.steepest_descent().max_steps(0));
		
		energy_match.check("total energy",        result.energy.total(),         -11.569230679378);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.513350460378);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.715154057693);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.974854033969);
		energy_match.check("external energy",     result.energy.external(),       -5.856397409990);
		energy_match.check("non-local energy",    result.energy.non_local(),       -0.593446829105);
		energy_match.check("XC energy",           result.energy.xc(),             -3.357693539830);
		energy_match.check("XC density integral", result.energy.nvxc(),           -3.698021814519);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),    -0.515173266198);
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}

	if(energy_match.fail()) return energy_match.fail();

	{
		auto result = ground_state::calculate(ions, electrons, options::theory{}.hartree_fock(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));

		energy_match.check("total energy",        result.energy.total(),          -9.788709725748);
		energy_match.check("kinetic energy",      result.energy.kinetic(),         8.151819376871);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),    -0.249826944617);
		energy_match.check("Hartree energy",      result.energy.hartree(),         0.880126063939);
		energy_match.check("external energy",     result.energy.external(),       -5.499063363813);
		energy_match.check("non-local energy",    result.energy.non_local(),       -0.510900262733);
		energy_match.check("XC energy",           result.energy.xc(),              0.000000000000);
		energy_match.check("XC density integral", result.energy.nvxc(),            0.000000000000);
		energy_match.check("HF exchange energy",  result.energy.exact_exchange(),    -2.075967411411);		
		energy_match.check("ion-ion energy",      result.energy.ion(),           -10.734724128603);
	}
	
	fftw_cleanup();
	
	return energy_match.fail();
	
}


