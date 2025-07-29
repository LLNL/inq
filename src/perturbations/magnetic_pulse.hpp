/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__MAGNETIC_PULSE
#define INQ__PERTURBATIONS__MAGNETIC_PULSE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Yifan Yao
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <magnitude/chirp.hpp>

namespace inq {
namespace perturbations {

class magnetic_pulse : public perturbations::none {

    vector3<double> magnetic_vector_;
    vector3<double> pulse_vector_;
    vector3<double> step_vector_;
    double frequency_;
    double width_;
    double chirp_;
    double time_shift_;
    double phase_shift_;

public:
    magnetic_pulse(vector3<quantity<magnitude::magfield>> magnetic_vector,
        vector3<quantity<magnitude::magfield>> pulse_vector,
        vector3<quantity<magnitude::magfield>> step_vector,
        quantity<magnitude::energy> frequency,
        quantity<magnitude::time> width,
        quantity<magnitude::chirp> chirp,
        quantity<magnitude::time> time_shift,
        double phase_shift):
        magnetic_vector_(vector3<double> {magnetic_vector[0].in_atomic_units(), magnetic_vector[1].in_atomic_units(), magnetic_vector[2].in_atomic_units()}),
        pulse_vector_(vector3<double> {pulse_vector[0].in_atomic_units(), pulse_vector[1].in_atomic_units(), pulse_vector[2].in_atomic_units()}),
        step_vector_(vector3<double> {step_vector[0].in_atomic_units(), step_vector[1].in_atomic_units(), step_vector[2].in_atomic_units()}),
        frequency_(frequency.in_atomic_units()),
        width_(width.in_atomic_units()),
        chirp_(chirp.in_atomic_units()),
        time_shift_(time_shift.in_atomic_units()),
        phase_shift_(phase_shift)
    {
    }

    auto has_magnetic_field() const {
        return true;
    }

    auto uniform_magnetic_field(double time) const {
        if (time > time_shift_) {
            return pulse_vector_*( cos(pow(time-time_shift_,2.0)*chirp_ + (time-time_shift_)*frequency_ + phase_shift_) * exp(-0.5*pow(time-time_shift_,2.0)/pow(width_,2.0)) - 1.0) + step_vector_ + magnetic_vector_;
        }
        else {
            return magnetic_vector_;
        }
    }

    template<typename MagneticField>
    void magnetic_field(const double time, MagneticField & magnetic) const {
        gpu::run(magnetic.basis().local_size(),
            [magnetic_ = begin(magnetic.linear()), mv = magnetic_vector_, pv = pulse_vector_, sv = step_vector_, t_ = time, t0_ = time_shift_, f_ = frequency_, phs_ = phase_shift_, ch_ = chirp_, w_ = width_] GPU_LAMBDA (auto ip){
                if (t_ > t0_) {
                    magnetic_[ip] += pv*( cos(pow(t_-t0_,2.0)*ch_ + (t_-t0_)*f_ + phs_) * exp(-0.5*pow(t_-t0_,2.0)/pow(w_,2.0)) - 1.0) + sv + mv;
                }
                else {
                    magnetic_[ip] += mv;
                }
            });
    }

};



}
}
#endif

#ifdef INQ_PERTURBATIONS_MAGNETIC_PULSE_UNIT_TEST
#undef INQ_PERTURBATIONS_MAGNETIC_PULSE_UNIT_TEST

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){
    
    parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

    perturbations::magnetic_pulse uniform_mag_pulse{{0.0_amu, 0.0_amu, 1.0_amu}, {0.0_amu, 0.0_amu, 0.0_amu}, {0.0_amu, 0.0_amu, 0.0_amu}, 1.0_ev, 1.0_attosecond, 1.0_invfs2, 100.0_attosecond, 0.0};

    CHECK(uniform_mag_pulse.has_magnetic_field());

}
#endif