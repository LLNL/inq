#include <inq/inq.hpp>
#include <stdio.h>

using namespace inq;
using namespace inq::magnitude;

inq::utils::match match(3.0e-4);

template <class EnvType>
auto H_energy(EnvType const & env) {
    // geometry
    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});

    // collinear
    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());

    ground_state::initial_guess(ions, electrons);

    auto result_col = ground_state::calculate(ions, electrons, options::theory{}.non_interacting(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
    result_col = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha).max_steps(0));

    // non collinear calculation
    auto electrons_2 = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());

    ground_state::initial_guess(ions, electrons_2);

    auto result_noncol = ground_state::calculate(ions, electrons_2, options::theory{}.non_interacting(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha));
    result_noncol = ground_state::calculate(ions, electrons_2, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1e-8_Ha).max_steps(0));

    match.check("total energy",           result_col.energy.total(),       result_noncol.energy.total());
    match.check("kinetic energy",         result_col.energy.kinetic(),     result_noncol.energy.kinetic());
    match.check("eigenvalues",            result_col.energy.eigenvalues(), result_noncol.energy.eigenvalues());
    match.check("Hartree energy",         result_col.energy.hartree(),     result_noncol.energy.hartree());
    match.check("External energy",        result_col.energy.external(),    result_noncol.energy.external());
    match.check("Non-local energy",       result_col.energy.non_local(),   result_noncol.energy.non_local());
    match.check("XC energy",              result_col.energy.xc(),          result_noncol.energy.xc());
    match.check("nvXC energy",            result_col.energy.nvxc(),        result_noncol.energy.nvxc());

}

int main (int argc, char ** argv) {
    auto env = input::environment{};

    // compute energy
    H_energy(env);
}