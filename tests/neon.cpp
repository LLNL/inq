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

	auto comm = input::environment::global().comm();

	utils::match energy_match(2.0e-5);

	interface::clear();
	
	//inq cell cubic 15.0 bohr finite
	interface::cell.cubic(15.0_b, 0);

	//inq cell
	interface::cell();

	//inq ions add Ne 0.0 0.0 0.0 bohr
	interface::ions.insert("Ne", {0.0_b, 0.0_b, 0.0_b});

  //inq ions
	interface::ions();
	
	//inq electrons extra_states 3
	interface::electrons.extra_states(3);
	
	//inq electrons cutoff 30.0 Ha
	interface::electrons.cutoff(30.0_Ha);

	//inq theory non_interacting
	interface::theory.non_interacting();

	//REAL SPACE PSEUDO
	//inq run ground_state
	interface::run.ground_state();
	
	energy_match &= interface::util.match(interface::result.energy_total()      , -61.861056649453, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_kinetic()    ,  35.765610219604, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_eigenvalues(), -61.861045337100, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_external()   , -79.509954154661, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_non_local()  , -18.116701402044, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_hartree()    ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_xc()         ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_nvxc()       ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_nvxc()       ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_ion()        ,   0.000000000000, 2.0e-5);
	
	//FOURIER SPACE PSEUDO
	//inq electrons fourier_pseudo
	interface::electrons.fourier_pseudo();

	//inq run ground_state
	interface::run.ground_state();
	
	energy_match &= interface::util.match(interface::result.energy_total()      , -61.861056649453, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_kinetic()    ,  35.765555684056, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_eigenvalues(), -61.861056649453, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_external()   , -79.509918897873, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_non_local()  , -18.116693435635, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_hartree()    ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_xc()         ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_nvxc()       ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_nvxc()       ,   0.000000000000, 2.0e-5);
	energy_match &= interface::util.match(interface::result.energy_ion()        ,   0.000000000000, 2.0e-5);
	
	return energy_match.fail();
}
