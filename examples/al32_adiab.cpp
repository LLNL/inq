/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa
*/

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <operations/io.hpp>
#include <utils/match.hpp>
#include <ground_state/initialize.hpp>
#include <ground_state/calculate.hpp>
#include <input/environment.hpp>
#include <input/parse_xyz.hpp>
#include <config/path.hpp>

#include <real_time/propagate.hpp>

#include<fstream>
#include <iomanip>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env(argc, argv);

	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
	
	utils::match energy_match(4.0e-6);

	std::vector<input::atom> geo = input::parse_xyz(config::path::unit_tests_data() + "al32.xyz");

	auto const L = 2*7.6524459;
	int const N = 200;
	double const dx = L/N;

	std::ofstream ofs{"al32_adiab.dat"}; 
	ofs<<"# distance (au), energy (au)\n"<< std::setprecision(10);

	input::config conf;
	
	conf.excess_charge = -1;
	conf.extra_states = 8;
	conf.temperature = 300.0_K;

	geo.emplace_back("H" | math::vector3<double>(0.00000, 1.91325, 1.91325));

	systems::ions ions(input::cell::cubic(L*1._b), geo);
	
	systems::electrons electrons(comm_world, ions, input::basis::cutoff_energy(25.0_Ha), conf);
	
	ground_state::initialize(ions, electrons);
	
	for(int n = 0; n != N; ++n){

		double const x = dx*n;

		ions.geo().coordinates()[ions.geo().num_atoms() - 1] = {0.00000 + x, 1.91325, 1.91325};
	
		auto result = ground_state::calculate(
			ions, electrons, input::interaction::pbe(), 
			inq::input::scf::steepest_descent() | inq::input::scf::scf_steps(400)
		);

		ofs<< x <<'\t'<< result.energy.total() << std::endl;
	}

	fftw_cleanup(); //required for valgrid
	
	return 0;
}

