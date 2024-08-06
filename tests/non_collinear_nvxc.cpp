#include <inq/inq.hpp>
#include <operations/integral.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <stdio.h>

using namespace inq;
using namespace inq::magnitude;

inq::utils::match match(3.0e-4);

template<class occupations_array_type, class field_set_type, typename VxcType, typename SpinDensityType>
auto compute_psi_vxc_psi(occupations_array_type const & occupations, field_set_type const & phi, VxcType const & vxc, SpinDensityType const & density) {
    
    int nst = phi.local_spinor_set_size();
    int npt = phi.basis().local_size();
    states::orbital_set<basis::real_space, double> vphi(phi.skeleton());

    inq::basis::field<basis::real_space, double> rfield(vphi.basis());

    double nvxc = 0.0;

    if (phi.spinors()){
        assert(std::get<1>(sizes(phi.spinor_array())) == phi.spinor_dim());
        assert(std::get<2>(sizes(phi.spinor_array())) == phi.local_spinor_set_size());
        //v = begin(vxc.matrix())
        //if (vxc.set_size() == 2) {
        //    gpu::run(phi.local_spinor_set_size(), phi.basis().local_size(),
        //            [ph = begin(phi.spinor_array()), vph = begin(vphi.matrix())] GPU_LAMBDA (auto ist, auto ip){
        //                vph[ip][ist] += norm(ph[ip][0][ist]);
        //                vph[ip][ist] += norm(ph[ip][1][ist]);
        //            });
        //    for (int ist=0; ist < nst; ist++) {
        //        basis::field_set<basis::real_space, double> ph(density.skeleton());
        //        ph.fill(0.0);
        //        gpu::run(phi.basis().local_size(),
        //            [p = begin(ph.matrix()), vph = begin(vphi.matrix()), ist] GPU_LAMBDA (auto ip){
        //                p[ip] += vph[ip][ist];
        //            });
                //auto r = operations::integral(ph);
                //nvxc += occupations[ist] * r;
        //    }
        //}

        //gpu::run(phi.local_spinor_set_size(), phi.basis().local_size(),
        //            [ph = begin(phi.spinor_array()), vph = begin(vphi.matrix())] GPU_LAMBDA (auto ist, auto ip){
        //                vph[ip][ist] += norm(ph[ip][0][ist]);
        //                vph[ip][ist] += norm(ph[ip][1][ist]);
        //            });
    }

    //nvxc += operations::integral_sum(vphi);
    
    return nvxc;
}

template <class EnvType>
auto compare_nvxc(EnvType const & env, systems::ions const & ions) {

    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
    ground_state::initial_guess(ions, electrons);
    auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(100));
    auto nvxc = result.energy.nvxc();
    std::cout << " NVXC :  " << nvxc << std::endl;

    auto core_density_ = electrons.atomic_pot().nlcc_density(electrons.states_comm(), electrons.spin_density().basis(), ions);
    basis::field_set<basis::real_space, double> vxc(electrons.spin_density().skeleton());
    hamiltonian::xc_term xc_(options::theory{}.lda(), electrons.spin_density().set_size());
    xc_.compute_vxc(electrons.spin_density(), core_density_, vxc);

    int iphi = 0;
    auto nvxc_2 = 0.0;
    for (auto & phi : electrons.kpin()) {
        nvxc_2 += compute_psi_vxc_psi(electrons.occupations()[iphi], phi, vxc, electrons.spin_density());
        iphi++;
    }
    std::cout << iphi << std::endl;
    std::cout << " NVXC2 : " << nvxc_2 << std::endl;
}


int main (int argc, char ** argv) {
    auto env = input::environment{};

    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
    compare_nvxc(env, ions);
}