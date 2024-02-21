/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__RAMPLASER
#define INQ__PERTURBATIONS__RAMPLASER

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Yifan Yao
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

class ramplaser : public perturbations::none {

	vector3<double, cartesian> polarization_;
	double frequency_;
	double rampstart_;
	double rampwidth_;
	gauge gauge_;

public:
	
	ramplaser(vector3<double, cartesian> polarization, quantity<magnitude::energy> frequency, quantity<magnitude::time> rampstart, quantity<magnitude::time> rampwidth, gauge arg_gauge = gauge::automatic):
		polarization_(polarization),
		frequency_(frequency.in_atomic_units()),
		rampstart_(rampstart.in_atomic_units()),
		rampwidth_(rampwidth.in_atomic_units()),
		gauge_(arg_gauge) {
		if(gauge_ == gauge::automatic) gauge_ = gauge::velocity;
	}
	
	auto has_uniform_electric_field() const {
		return gauge_ == gauge::length;
	}
	
	auto uniform_electric_field(double time) const {
		double coshfactor = cosh(rampwidth_*(time-rampstart_));
		return polarization_*sin(time*frequency_) *0.5*(tanh((time-rampstart_)/rampwidth_)+1.0) - polarization_/frequency_*(cos(time*frequency_) - 1.0) * 0.5/rampwidth_/coshfactor/coshfactor;
	}
	
	auto has_uniform_vector_potential() const {
		return gauge_ == gauge::velocity;
	}
	
	auto uniform_vector_potential(double time) const {
		//E=-1/c*dA/dt
		return polarization_/frequency_*(cos(time*frequency_) - 1.0) *0.5*(tanh((time-rampstart_)/rampwidth_)+1.0);
	}

	template<class OStream>
	friend OStream & operator<<(OStream & out, ramplaser const & self){
		using namespace magnitude;

		auto freq_ev = self.frequency_/in_atomic_units(1.0_eV);
		out << "Ramp-laser:\n";
		out << "  polarization [a.u.] = " << self.polarization_ << "\n";
		out << "  frequency           = " << self.frequency_ << " Ha | " << freq_ev << " eV | " << freq_ev*241.7991 << " THz | " << 1239.84193/freq_ev << " nm" << std::endl;
		out << "  gauge               = " << self.gauge_ << "\n";
		return out;
	}
	
};



}
}
#endif

#ifdef INQ_PERTURBATIONS_RAMPLASER_UNIT_TEST
#undef INQ_PERTURBATIONS_RAMPLASER_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	perturbations::ramplaser rlas({1.0, 0.0, 0.0}, 1.0_eV, 0.0_fs, 1.0_fs);

	std::cout << rlas;
	
	CHECK(rlas.has_uniform_electric_field());

	SECTION("ramplaser velocity gauge"){
		perturbations::ramplaser ramp_vector_potential({0.1, 0.0, 0.0}, 1.0_eV, 0.0_fs, 1.0_fs, perturbations::gauge::velocity);
		CHECK(ramp_vector_potential.has_uniform_vector_potential());
		CHECK(not ramp_vector_potential.has_uniform_electric_field());
		CHECK(ramp_vector_potential.uniform_vector_potential(0.0)[0] == 0.0);
		CHECK(ramp_vector_potential.uniform_vector_potential(0.0)[2] == 0.0);
	}

	SECTION("ramplaser length gauge"){
		perturbations::ramplaser ramp_E_field({1.0, 0.0, 0.0}, 1.0_eV, 0.0_fs, 1.0_fs, perturbations::gauge::length);
		CHECK(ramp_E_field.has_uniform_electric_field());
		CHECK(not ramp_E_field.has_uniform_vector_potential());
		CHECK(ramp_E_field.uniform_electric_field(0.0)[0] == 0.0);
	}
}
#endif
