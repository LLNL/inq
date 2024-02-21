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
	
	laser(vector3<double, cartesian> polarization, quantity<magnitude::energy> frequency,gauge arg_gauge = gauge::length):
		polarization_(polarization),
		frequency_(frequency.in_atomic_units()),
		gauge_(arg_gauge) {
		assert(gauge_ != gauge::mixed);
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

	template<class OStream>
	friend OStream & operator<<(OStream & out, laser const & self){
		using namespace magnitude;

		auto freq_ev = self.frequency_*27.211383;
		
		out << "Frequency :    " << self.frequency_ << " Ha" << std::endl;
		out << "               " << freq_ev << " eV" << std::endl;
		out << "               " << freq_ev*241.7991 << " THz" << std::endl;
		out << "               " << 1239.84193/freq_ev << " nm" << std::endl;

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

	perturbations::laser las({1.0, 0.0, 0.0}, 1.0_eV);

	std::cout << las;
	
	CHECK(las.has_uniform_electric_field());

	SECTION("velocity gauge"){
		perturbations::laser vector_potential({0.1, 0.0, 0.0}, 1.0_eV,perturbations::gauge::velocity);
		CHECK(vector_potential.has_uniform_vector_potential());
		CHECK(not vector_potential.has_uniform_electric_field());
		CHECK(vector_potential.uniform_vector_potential(0.0)[0] == 0.0);
		CHECK(vector_potential.uniform_vector_potential(0.0)[2] == 0.0);
	}

	SECTION("length gauge"){
		perturbations::laser E_field({1.0, 0.0, 0.0}, 1.0_eV);
		CHECK(E_field.has_uniform_electric_field());
		CHECK(not E_field.has_uniform_vector_potential());
		CHECK(E_field.uniform_electric_field(0.0)[0] == 0.0);
	}

}
#endif
