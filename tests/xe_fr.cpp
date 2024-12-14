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

	systems::ions ions(systems::cell::cubic(10.0_b).finite());
	ions.species_list().pseudopotentials() = pseudo::set_id::pseudodojo_rel_pbe();
	//	ions.species_list().insert(ionic::species("Xe").pseudo_file(config::path::unit_tests_data() + "Xe_fr.UPF.gz"));
	ions.insert("Xe", {0.0_b, 0.0_b, 0.0_b});
	
	systems::electrons electrons(ions, options::electrons{}.cutoff(40.0_Ha).spin_non_collinear().extra_states(3));
	ground_state::initial_guess(ions, electrons);

	//we use LDA for now since PBE has some convergence issues
	auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), options::ground_state{}.energy_tolerance(1e-8_Ha));

	auto all_eigenvalues = parallel::gather(+electrons.eigenvalues().flatted(), electrons.kpin_states_comm(), 0);
	
	if(electrons.kpin_states_comm().root()){
		energy_match.check("eigenvalue  0",    all_eigenvalues[ 0],  -0.712820889809);
		energy_match.check("eigenvalue  1",    all_eigenvalues[ 1],  -0.712820889809);
		energy_match.check("eigenvalue  2",    all_eigenvalues[ 2],  -0.322520169668);
		energy_match.check("eigenvalue  3",    all_eigenvalues[ 3],  -0.322520169668);
		energy_match.check("eigenvalue  4",    all_eigenvalues[ 4],  -0.277257682172);
		energy_match.check("eigenvalue  5",    all_eigenvalues[ 5],  -0.277257682172);
		energy_match.check("eigenvalue  6",    all_eigenvalues[ 6],  -0.277257669711);
		energy_match.check("eigenvalue  7",    all_eigenvalues[ 7],  -0.277257669711);
		energy_match.check("eigenvalue  8",    all_eigenvalues[ 8],  -0.063230692730);
		energy_match.check("eigenvalue  9",    all_eigenvalues[ 9],  -0.063230692417);
		energy_match.check("eigenvalue 10",    all_eigenvalues[10],   0.046161112021);
	}
	
	energy_match.check("total energy",        result.energy.total(),          -18.661808075647);
	energy_match.check("kinetic energy",      result.energy.kinetic(),          4.169408422510);
	energy_match.check("eigenvalues",         result.energy.eigenvalues(),     -3.179712822721);
	energy_match.check("Hartree energy",      result.energy.hartree(),         13.167234189885);
	energy_match.check("external energy",     result.energy.external(),       -32.516230615704);
	energy_match.check("non-local energy",    result.energy.non_local(),        2.758847694594);
	energy_match.check("XC energy",           result.energy.xc(),              -6.241067766933);
	energy_match.check("XC density integral", result.energy.nvxc(),            -3.926206703892);
	energy_match.check("HF exchange energy",  result.energy.exact_exchange(),   0.000000000000);
	energy_match.check("ion-ion energy",      result.energy.ion(),              0.000000000000);
	
	return energy_match.fail();
	
}
