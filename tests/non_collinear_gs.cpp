#include <inq/inq.hpp>
#include <stdio.h>

using namespace inq;
using namespace inq::magnitude;

inq::utils::match match(3.0e-4);

template <class EnvType>
auto compute_GS(EnvType const & env, systems::ions const & ions) {

    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
    ground_state::initial_guess(ions, electrons);
    auto result_col = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(100));
    auto mag = observables::total_magnetization(electrons.spin_density());

    auto electrons_nc = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
    ground_state::initial_guess(ions, electrons_nc);
    auto result_nc = ground_state::calculate(ions, electrons_nc, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(100));
    auto mag_nc = observables::total_magnetization(electrons_nc.spin_density());

    match.check("total magnetization",        mag.length(),      mag_nc.length());
}


int main (int argc, char ** argv) {
    auto env = input::environment{};

    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
    compute_GS(env, ions);

    auto d = 1.21_A;
    ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("O", {0.0_b, 0.0_b, d/2});
    ions.insert("O", {0.0_b, 0.0_b,-d/2});
    compute_GS(env, ions);
}