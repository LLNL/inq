/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <fftw3.h>

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq::input;
	using namespace inq::systems;
	using namespace inq::magnitude;
	using inq::vector3;
	
	inq::input::environment env(argc, argv);
	
	inq::utils::match match(3.0e-4);
	
	auto ions = inq::systems::ions::parse(inq::config::path::unit_tests_data() + "water.xyz", inq::ions::unit_cell::orthorhombic(12.0_b, 11.0_b, 10.0_b).finite());
	
	auto comm = boost::mpi3::environment::get_world_instance();
	auto parstates = comm.size();
	if(comm.size() == 3 or comm.size() == 5) parstates = 1;
	
	inq::systems::electrons electrons(env.par().states(parstates), ions, inq::input::config::cutoff(30.0_Ha));

	inq::ground_state::initial_guess(ions, electrons);

	auto scf_options = scf::energy_tolerance(1.0e-9_Ha) | scf::broyden_mixing();
	auto result = inq::ground_state::calculate(ions, electrons, interaction::lda(), scf_options);
	
	match.check("total energy",        result.energy.total(),       -17.604152928274);
	match.check("kinetic energy",      result.energy.kinetic(),      12.055671278976);
	match.check("eigenvalues",         result.energy.eigenvalues(),  -4.066524514529);
	match.check("Hartree energy",      result.energy.hartree(),      21.255096093096);
	match.check("external energy",     result.energy.external(),    -50.405054484947);
	match.check("non-local energy",    result.energy.nonlocal(),     -2.732683697594);
	match.check("XC energy",           result.energy.xc(),           -4.762305356613);
	match.check("XC density integral", result.energy.nvxc(),         -5.494649797157);
	match.check("HF exchange energy",  result.energy.hf_exchange(),   0.000000000000);
	match.check("ion-ion energy",      result.energy.ion(),           6.985123238808);

	std::cout << result.dipole << std::endl;
	
	match.check("dipole x", result.dipole[0], -0.000304523);
	match.check("dipole y", result.dipole[1], -0.724304);
	match.check("dipole z", result.dipole[2], -2.78695e-05);
	
	electrons.save("h2o_restart");
	
	fftw_cleanup(); //required for valgrid
	
	return match.fail();

}

