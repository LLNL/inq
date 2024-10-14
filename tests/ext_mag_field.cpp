/* -*- indent-tabs-mode: t -*- */

#include <inq/inq.hpp>
#include <perturbations/magnetic_field.hpp>

using namespace inq;
using namespace inq::magnitude;

int main (int argc, char ** argv){
    auto env = input::environment{};

    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});

    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
    ground_state::initial_guess(ions, electrons);
    perturbations::magnetic_field bfield{{0.0, 0.0, 1.0}};
    auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), bfield);
}