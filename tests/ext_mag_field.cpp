/* -*- indent-tabs-mode: t -*- */

#include <inq/inq.hpp>
#include <perturbations/magnetic.hpp>

using namespace inq;
using namespace inq::magnitude;

inq::utils::match match(3.0e-4);

int main (int argc, char ** argv){
    auto env = input::environment{};

    auto d = 1.21_A;
    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("O", {0.0_b, 0.0_b, d/2});
    ions.insert("O", {0.0_b, 0.0_b,-d/2});

    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(10).spin_polarized());
    ground_state::initial_guess(ions, electrons);
    perturbations::magnetic B{{0.0, 0.0, 0.4}};
    auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), B);
    auto mag = observables::total_magnetization(electrons.spin_density());
    std::cout << mag << std::endl;
    match.check("magnetization direction",        mag/mag.length(), {0.0, 0.0, 1.0});

    perturbations::magnetic B2{{0.0, 0.0, -0.4}};
    result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), B2);
    auto mag2 = observables::total_magnetization(electrons.spin_density());
    std::cout << mag2 << std::endl;
    match.check("magnetization direction",        mag2/mag2.length(), {0.0, 0.0, -1.0});

    perturbations::magnetic B3{{0.001, 0.0, -0.4}};
    electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(10).spin_non_collinear());
    ground_state::initial_guess(ions, electrons);
    result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.01), B3);
    auto mag3 = observables::total_magnetization(electrons.spin_density());
    std::cout << mag3 << std::endl;
    vector3 e_v = {0.001, 0.0, -0.4};
    e_v = e_v/sqrt(norm(e_v));
    match.check("magnetization direction",        mag3/mag3.length(), e_v);

    return match.fail();
}