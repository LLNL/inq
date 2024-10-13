#include <inq/inq.hpp>
#include <operations/integral.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <stdio.h>

using namespace inq;
using namespace inq::magnitude;

inq::utils::match match(3.0e-4);

template<class occupations_array_type, class field_set_type> 
double compute_psi_psi(occupations_array_type const & occupations, field_set_type const & phi) {

    basis::field<basis::real_space, double> rfield(phi.basis());
    rfield.fill(0.0);
    assert(std::get<1>(sizes(phi.spinor_array())) == phi.spinor_dim());
    assert(std::get<2>(sizes(phi.spinor_array())) == phi.local_spinor_set_size());

    if (phi.spinor_dim() == 1) {
        gpu::run(phi.local_spinor_set_size(), phi.basis().local_size(),
            [ph = begin(phi.matrix()), rf = begin(rfield.linear()), occ = begin(occupations)] GPU_LAMBDA (auto ist, auto ip){
                rf[ip] += occ[ist] * norm(ph[ip][ist]);
            });
    }
    else {
        gpu::run(phi.local_spinor_set_size(), phi.basis().local_size(),
            [ph = begin(phi.spinor_array()), rf = begin(rfield.linear()), occ = begin(occupations)] GPU_LAMBDA (auto ist, auto ip){
                rf[ip] += occ[ist] * norm(ph[ip][0][ist]);
                rf[ip] += occ[ist] * norm(ph[ip][1][ist]);
            });
    }
    auto totn = operations::integral(rfield);
    return totn;
}

template <class EnvType>
auto compare_nvxc(EnvType const & env, systems::ions const & ions, systems::electrons & electrons) {

    ground_state::initial_guess(ions, electrons);
    auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000));
    auto nvxc = result.energy.nvxc();
    auto exc = result.energy.xc();

    auto core_density_ = electrons.atomic_pot().nlcc_density(electrons.states_comm(), electrons.spin_density().basis(), ions);
    hamiltonian::xc_term xc_(options::theory{}.lda(), electrons.spin_density().set_size());

    auto nvxc_2 = 0.0;
    auto exc_2 = 0.0;
    xc_.eval_psi_vxc_psi(electrons.kpin_states_comm(), core_density_, electrons.spin_density(), electrons.occupations(), electrons.kpin(), nvxc_2, exc_2);

    auto totn = 0.0;
    int iphi = 0;
    for (auto & phi : electrons.kpin()) {
        totn += compute_psi_psi(electrons.occupations()[iphi], phi);
        iphi++;
    }

    match.check("electrons number",      electrons.states().num_electrons(),       totn);
    match.check("nvxc",                  nvxc_2,                                   nvxc);
    match.check("xc",                    exc_2,                                    exc);
}


int main (int argc, char ** argv) {
    auto env = input::environment{};

    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_unpolarized());
    compare_nvxc(env, ions, electrons);

    auto electrons_2 = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
    compare_nvxc(env, ions, electrons_2);

    auto electrons_3 = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
    compare_nvxc(env, ions, electrons_3);
}
