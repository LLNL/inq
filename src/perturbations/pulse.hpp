/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__PULSE
#define INQ__PERTURBATIONS__PULSE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Yifan Yao
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>
#include <magnitude/chirp.hpp>
#include <perturbations/gauge.hpp>
#include <perturbations/none.hpp>

#include <cmath>

namespace inq {
namespace perturbations {

class pulse : public perturbations::none {

	vector3<double, cartesian> polarization_;
	double frequency_;
	double width_;
	double chirp_;
	double time_shift_;
	double phase_shift_;
	
public:
	pulse(vector3<double, cartesian> polarization,
            quantity<magnitude::energy> frequency, 
            quantity<magnitude::time>  width,  
            quantity<magnitude::chirp> chirp, 
            quantity<magnitude::time> time_shift, 
            double phase_shift,
            gauge arg_gauge = gauge::length):
		polarization_(polarization),
		frequency_(frequency.in_atomic_units()),
                width_(width.in_atomic_units()),
                chirp_(chirp.in_atomic_units()),
                time_shift_(time_shift.in_atomic_units()),
                phase_shift_(phase_shift),
		gauge_(arg_gauge) {
		if(gauge_ == gauge::automatic) gauge_ = gauge::velocity;
	}
	
	auto has_uniform_electric_field() const {
		return gauge_ == gauge::length;
	}
	
	auto uniform_electric_field(double time) const {
		//E=-1/c*dA/dt
          double phase  = pow(time-time_shift_,2.0)*chirp_ + (time-time_shift_)*frequency_ + phase_shift_;
          double exponent = -0.5*pow(time-time_shift_,2.0)/pow(width_,2.0);
	  return polarization_/frequency_ * (2.0*(time-time_shift_)*chirp_+frequency_) * sin(phase) * exp(exponent) 
                  - polarization_/frequency_* cos(phase) * exp(exponent) * (-1.0*(time-time_shift_)/pow(width_,2.0)) ;
	}
	
	auto has_uniform_vector_potential() const {
		return gauge_ == gauge::velocity;
	}
	
	auto uniform_vector_potential(double time) const {
		//E=-1/c*dA/dt
		return polarization_/frequency_*( cos(pow(time-time_shift_,2.0)*chirp_ + (time-time_shift_)*frequency_ + phase_shift_) * exp(-0.5*pow(time-time_shift_,2.0)/pow(width_,2.0)) - 1.0);
	}

	template <typename OutputStream>
	void print_info(OutputStream & out) {
		auto freq_ev = frequency_*27.211383;
		
		out << "Frequency :    " << frequency_ << " Ha" << std::endl;
		out << "               " << freq_ev << " eV" << std::endl;
		out << "               " << freq_ev*241.7991 << " THz" << std::endl;
		out << "               " << 1239.84193/freq_ev << " nm" << std::endl;
		
	}
	
private:
	gauge gauge_;

};

}
}
#endif

#ifdef INQ_PERTURBATIONS_PULSE_UNIT_TEST
#undef INQ_PERTURBATIONS_PULSE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

        perturbations::pulse las({1.0, 0.0, 0.0}, 1.0_eV, 1.0_fs, 0.0_invfs2, 0.0_fs, M_PI, perturbations::gauge::length);

	SECTION("ramplaser velocity gauge"){
                perturbations::pulse las({1.0, 0.0, 0.0}, 1.0_eV, 1.0_fs, 0.0_invfs2, 0.0_fs, 0.0, perturbations::gauge::velocity);
		CHECK(las.has_uniform_vector_potential());
		CHECK(not las.has_uniform_electric_field());
		CHECK(las.uniform_vector_potential(0.0)[0] == 0.0);
		CHECK(las.uniform_vector_potential(0.0)[2] == 0.0);
	}

}
#endif
