/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <utils/match.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>

#include <input/environment.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
			
	input::environment env(argc, argv);
	
	utils::match energy_match(1.0e-5);

	{
		systems::ions ions(ions::unit_cell::cubic(10.0_b));
		systems::electrons electrons(env.par(), ions, input::kpoints::grid({1, 1, 3}), input::config{}.cutoff(30.0_Ha).temperature(300.0_K).extra_states(2).excess_charge(14.0));
		
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::lda(), inq::input::scf::energy_tolerance(1e-9_Ha));
		
		energy_match.check("total energy",        result.energy.total()      , -0.567967321401);
		energy_match.check("kinetic energy",      result.energy.kinetic()    ,  2.485678165423);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(), -1.488505428934);
		energy_match.check("hartree",             result.energy.hartree()    ,  0.000000732036);	
		energy_match.check("XC energy",           result.energy.xc()         , -3.053646218860);
		energy_match.check("XC density integral", result.energy.nvxc()       , -3.974185058430);
	}

	{
		auto a = 10.0_b;
		systems::ions ions(ions::unit_cell::lattice({a/sqrt(2.0), a/2.0, a/2.0}, {-a/sqrt(2), a/2.0, a/2.0}, {0.0_b, -a/sqrt(2.0), a/sqrt(2.0)}));
		systems::electrons electrons(env.par(), ions, input::config{}.cutoff(30.0_Ha).temperature(300.0_K).extra_states(2).excess_charge(14.0), input::kpoints::grid({1, 1, 3}));
		
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::lda(), inq::input::scf::energy_tolerance(1e-9_Ha));
		
		energy_match.check("total energy",        result.energy.total()      , -0.567967370592);
		energy_match.check("kinetic energy",      result.energy.kinetic()    ,  2.485678162550);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(), -1.488505616755);
		energy_match.check("hartree",             result.energy.hartree()    ,  0.000000551815);	
		energy_match.check("XC energy",           result.energy.xc()         , -3.053646084957);
		energy_match.check("XC density integral", result.energy.nvxc()       , -3.974184882934);
	}

	{
		auto a = 10.0_b;
		systems::ions ions(ions::unit_cell::lattice({0.0_b, a/2.0, a/2.0}, {a/2.0, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}));
		systems::electrons electrons(env.par(), ions, input::config{}.cutoff(30.0_Ha).temperature(300.0_K).extra_states(2).excess_charge(18.0), input::kpoints::grid({1, 1, 1}, false));
		
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons, input::interaction::lda(), inq::input::scf::energy_tolerance(1e-9_Ha));
		
		energy_match.check("total energy",        result.energy.total()      ,  3.023858102368);
		energy_match.check("kinetic energy",      result.energy.kinetic()    ,  9.474820227644);
		energy_match.check("eigenvalues",         result.energy.eigenvalues(),  1.054657729496);
		energy_match.check("hartree",             result.energy.hartree()    ,  0.000000000177);	
		energy_match.check("XC energy",           result.energy.xc()         , -6.450962125453);
		energy_match.check("XC density integral", result.energy.nvxc()       , -8.420162498501);
	}
		
	return energy_match.fail();
	
}
