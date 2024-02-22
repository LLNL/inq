/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__LASER
#define INQ__PERTURBATIONS__LASER

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Yifan Yao
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>
#include <perturbations/gauge.hpp>
#include <perturbations/none.hpp>

namespace inq {
namespace perturbations {

class laser : public perturbations::none {

	vector3<double, cartesian> polarization_;
	double frequency_;
	gauge gauge_;

public:
	
	laser(vector3<double, cartesian> polarization, quantity<magnitude::energy> frequency, gauge arg_gauge = gauge::automatic):
		polarization_(polarization),
		frequency_(frequency.in_atomic_units()),
		gauge_(arg_gauge) {

		if(gauge_ == gauge::automatic) gauge_ = gauge::velocity;
	}
	
	auto has_uniform_electric_field() const {
		return gauge_ == gauge::length;
	}
	
	auto uniform_electric_field(double time) const {
		return polarization_*sin(time*frequency_);
	}
	
	auto has_uniform_vector_potential() const {
		return gauge_ == gauge::velocity;
	}
	
	auto uniform_vector_potential(double time) const {
		//E=-1/c*dA/dt
		return polarization_/frequency_*(cos(time*frequency_) - 1.0);
	}

	void save(parallel::communicator & comm, std::string const & dirname) const {
		auto error_message = "INQ error: Cannot save the perturbations::laser to directory '" + dirname + "'.";
		
		utils::create_directory(comm, dirname);
		utils::save_value(comm, dirname + "/polarization",  polarization_,  error_message);
		utils::save_value(comm, dirname + "/frequency",     frequency_,     error_message);
		utils::save_value(comm, dirname + "/gauge",         gauge_,         error_message);
	}

	static auto load(std::string const & dirname) {
		using namespace magnitude;
		
    auto error_message = "INQ error: Cannot load perturbations::laser from directory '" + dirname + "'.";

		vector3<double> pol;
		double freq;
		gauge gau;
	
    utils::load_value(dirname + "/polarization",   pol,   error_message);
    utils::load_value(dirname + "/frequency",      freq,  error_message);
    utils::load_value(dirname + "/gauge",          gau,   error_message);
    
    return laser(pol, freq*1.0_Ha, gau);
	}

	template<class OStream>
	friend OStream & operator<<(OStream & out, laser const & self){
		using namespace magnitude;

		auto freq_ev = self.frequency_/in_atomic_units(1.0_eV);
		out << "Laser:\n";
		out << "  polarization [a.u.] = " << self.polarization_ << "\n";
		out << "  frequency           = " << self.frequency_ << " Ha | " << freq_ev << " eV | " << freq_ev*241.7991 << " THz | " << 1239.84193/freq_ev << " nm" << std::endl;
		out << "  gauge               = " << self.gauge_ << "\n";
		return out;
	}

};


}
}
#endif

#ifdef INQ_PERTURBATIONS_LASER_UNIT_TEST
#undef INQ_PERTURBATIONS_LASER_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
 
	SECTION("length gauge"){
		perturbations::laser las({1.0, 0.0, 0.0}, 1.0_eV, perturbations::gauge::length);

		std::cout << las;
	
		CHECK(las.has_uniform_electric_field());
		CHECK(not las.has_uniform_vector_potential());

		CHECK(las.uniform_electric_field(0.0)[0] == 0.0_a);
		CHECK(las.uniform_electric_field(0.0)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(0.0)[2] == 0.0_a);

		CHECK(las.uniform_electric_field(0.1)[0] == 0.0036749239_a);
		CHECK(las.uniform_electric_field(0.1)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(0.1)[2] == 0.0_a);

		CHECK(las.uniform_electric_field(0.5)[0] == 0.0183736271_a);
		CHECK(las.uniform_electric_field(0.5)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(0.5)[2] == 0.0_a);
		
		CHECK(las.uniform_electric_field(0.7)[0] == 0.0257216884_a);
		CHECK(las.uniform_electric_field(0.7)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(0.7)[2] == 0.0_a);
		
		CHECK(las.uniform_electric_field(1.4)[0] == 0.0514263564_a);
		CHECK(las.uniform_electric_field(1.4)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(1.4)[2] == 0.0_a);
		
		CHECK(las.uniform_electric_field(10.0)[0] == 0.3592771603_a);
		CHECK(las.uniform_electric_field(10.0)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(10.0)[2] == 0.0_a);
				
		CHECK(las.uniform_electric_field(70.0)[0] == 0.5389078998_a);
		CHECK(las.uniform_electric_field(70.0)[1] == 0.0_a);
		CHECK(las.uniform_electric_field(70.0)[2] == 0.0_a);

		las.save(comm, "save_laser_length");
		auto read_las = perturbations::laser::load("save_laser_length");

		CHECK(read_las.uniform_electric_field(0.0)[0] == 0.0_a);
		CHECK(read_las.uniform_electric_field(0.0)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(0.0)[2] == 0.0_a);

		CHECK(read_las.uniform_electric_field(0.1)[0] == 0.0036749239_a);
		CHECK(read_las.uniform_electric_field(0.1)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(0.1)[2] == 0.0_a);

		CHECK(read_las.uniform_electric_field(0.5)[0] == 0.0183736271_a);
		CHECK(read_las.uniform_electric_field(0.5)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(0.5)[2] == 0.0_a);
		
		CHECK(read_las.uniform_electric_field(0.7)[0] == 0.0257216884_a);
		CHECK(read_las.uniform_electric_field(0.7)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(0.7)[2] == 0.0_a);
		
		CHECK(read_las.uniform_electric_field(1.4)[0] == 0.0514263564_a);
		CHECK(read_las.uniform_electric_field(1.4)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(1.4)[2] == 0.0_a);
		
		CHECK(read_las.uniform_electric_field(10.0)[0] == 0.3592771603_a);
		CHECK(read_las.uniform_electric_field(10.0)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(10.0)[2] == 0.0_a);
				
		CHECK(read_las.uniform_electric_field(70.0)[0] == 0.5389078998_a);
		CHECK(read_las.uniform_electric_field(70.0)[1] == 0.0_a);
		CHECK(read_las.uniform_electric_field(70.0)[2] == 0.0_a);

	}

	SECTION("velocity gauge"){
		perturbations::laser las({0.0, -0.5, 0.5}, 1.0_eV, perturbations::gauge::velocity);

		CHECK(las.has_uniform_vector_potential());
		CHECK(not las.has_uniform_electric_field());

		CHECK(las.uniform_vector_potential(0.0)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(0.0)[1] ==  0.0_a);
		CHECK(las.uniform_vector_potential(0.0)[2] ==  0.0_a);

		CHECK(las.uniform_vector_potential(0.1)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(0.1)[1] ==  0.0000918732_a);
		CHECK(las.uniform_vector_potential(0.1)[2] == -0.0000918732_a);

		CHECK(las.uniform_vector_potential(0.5)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(0.5)[1] ==  0.002296768_a);
		CHECK(las.uniform_vector_potential(0.5)[2] == -0.002296768_a);
		
		CHECK(las.uniform_vector_potential(0.7)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(0.7)[1] ==  0.0045015437_a);
		CHECK(las.uniform_vector_potential(0.7)[2] == -0.0045015437_a);
		
		CHECK(las.uniform_vector_potential(1.4)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(1.4)[1] ==  0.0180031961_a);
		CHECK(las.uniform_vector_potential(1.4)[2] == -0.0180031961_a);
		
		CHECK(las.uniform_vector_potential(10.0)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(10.0)[1] ==  0.9084398165_a);
		CHECK(las.uniform_vector_potential(10.0)[2] == -0.9084398165_a);
				
		CHECK(las.uniform_vector_potential(70.0)[0] ==  0.0_a);
		CHECK(las.uniform_vector_potential(70.0)[1] ==  25.0666486301_a);
		CHECK(las.uniform_vector_potential(70.0)[2] == -25.0666486301_a);

		las.save(comm, "save_laser_velocity");
		auto read_las = perturbations::laser::load("save_laser_velocity");

		CHECK(read_las.uniform_vector_potential(0.0)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(0.0)[1] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(0.0)[2] ==  0.0_a);

		CHECK(read_las.uniform_vector_potential(0.1)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(0.1)[1] ==  0.0000918732_a);
		CHECK(read_las.uniform_vector_potential(0.1)[2] == -0.0000918732_a);

		CHECK(read_las.uniform_vector_potential(0.5)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(0.5)[1] ==  0.002296768_a);
		CHECK(read_las.uniform_vector_potential(0.5)[2] == -0.002296768_a);
		
		CHECK(read_las.uniform_vector_potential(0.7)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(0.7)[1] ==  0.0045015437_a);
		CHECK(read_las.uniform_vector_potential(0.7)[2] == -0.0045015437_a);
		
		CHECK(read_las.uniform_vector_potential(1.4)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(1.4)[1] ==  0.0180031961_a);
		CHECK(read_las.uniform_vector_potential(1.4)[2] == -0.0180031961_a);
		
		CHECK(read_las.uniform_vector_potential(10.0)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(10.0)[1] ==  0.9084398165_a);
		CHECK(read_las.uniform_vector_potential(10.0)[2] == -0.9084398165_a);
				
		CHECK(read_las.uniform_vector_potential(70.0)[0] ==  0.0_a);
		CHECK(read_las.uniform_vector_potential(70.0)[1] ==  25.0666486301_a);
		CHECK(read_las.uniform_vector_potential(70.0)[2] == -25.0666486301_a);
		
	}

}
#endif
