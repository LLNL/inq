/* -*- indent-tabs-mode: t -*- */

#ifndef ZEEMAN_COUPLING_HPP
#define ZEEMAN_COUPLING_HPP

#include <inq_config.h>

namespace inq {
namespace hamiltonian {

class zeeman_coupling {

private:

    int spin_components_;

public:

    zeeman_coupling(int const spin_components):
        spin_components_(spin_components)
    {
        assert(spin_components_ > 1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template<typename SpinDensityType, typename MagneticField, typename VKSType>
    void operator()(SpinDensityType const & spin_density, MagneticField const & B, VKSType & vks, double & zeeman_ener) const {

        basis::field_set<basis::real_space, double> vz(vks.skeleton());
        vz.fill(0.0);

        assert(vz.set_size() == spin_components_);

        compute_vz(B, vz);

        process_potential(vz, vks);

        zeeman_ener += compute_zeeman_energy(spin_density, vz);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template<typename MagneticField, typename VZType>
    void compute_vz(MagneticField const & B, VZType & vz) const {

        if (vz.set_size() == 4) {
            gpu::run(vz.basis().local_size(),
                [v = begin(vz.matrix()), b = begin(B.linear())] GPU_LAMBDA (auto ip) {
                    v[ip][0] +=-b[ip][2];
                    v[ip][1] += b[ip][2];
                    v[ip][2] +=-b[ip][0];
                    v[ip][3] +=-b[ip][1];
                });
        }
        else {
            assert(vz.set_size() == 2);
            gpu::run(vz.basis().local_size(),
                [v = begin(vz.matrix()), b = begin(B.linear())] GPU_LAMBDA (auto ip) {
                    v[ip][0] +=-b[ip][2];
                    v[ip][1] += b[ip][2];
                });
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template <typename SpinDensityType, typename VZType>
    double compute_zeeman_energy(SpinDensityType const & spin_density, VZType & vz) const {

        auto zeeman_ener_ = 0.0;
        if (spin_density.set_size() == 4) {
            gpu::run(spin_density.local_set_size(), spin_density.basis().local_size(),
                [v = begin(vz.matrix())] GPU_LAMBDA (auto is, auto ip) {
                    if (is >= 2) v[ip][is] = 2.0*v[ip][is];
                });
        }
        zeeman_ener_ += operations::integral_product_sum(spin_density, vz);
        return zeeman_ener_;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template<typename VZType, typename VKSType>
    void process_potential(VZType const & vz, VKSType & vks) const {

        gpu::run(vz.local_set_size(), vz.basis().local_size(),
            [v = begin(vz.matrix()), vk = begin(vks.matrix())] GPU_LAMBDA (auto is, auto ip) {
                vk[ip][is] += v[ip][is];
            });
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template<class occupations_array_type, class field_set_type, typename VZType, typename RFType>
    void compute_psi_vz_psi_ofr(occupations_array_type const & occupations, field_set_type const & phi, VZType const & vz, RFType & rfield) {

        assert(std::get<1>(sizes(phi.spinor_array())) == phi.spinor_dim());
        assert(std::get<2>(sizes(phi.spinor_array())) == phi.local_spinor_set_size());

        if (vz.set_size() == 2){
            gpu::run(phi.local_set_size(), phi.basis().local_size(),
                [ph = begin(phi.matrix()), rf = begin(rfield.linear()), v = begin(vz.matrix()), occ = begin(occupations), spi = phi.spin_index()] GPU_LAMBDA (auto ist, auto ip) {
                    rf[ip] += occ[ist]*v[ip][spi]*norm(ph[ip][ist]);
                });
        }
        else {
            assert(vz.set_size() == 4);
            gpu::run(phi.local_spinor_set_size(), phi.basis().local_size(),
                [ph = begin(phi.spinor_array()), rf = begin(rfield.linear()), v = begin(vz.matrix()), occ = begin(occupations)] GPU_LAMBDA (auto ist, auto ip) {
                    auto offdiag = v[ip][2] + complex{0.0, 1.0}*v[ip][3];
                    auto cross = 2.0*occ[ist]*real(offdiag*ph[ip][1][ist]*conj(ph[ip][0][ist]));
                    rf[ip] += occ[ist]*v[ip][0]*norm(ph[ip][0][ist]);
                    rf[ip] += occ[ist]*v[ip][1]*norm(ph[ip][1][ist]);
                    rf[ip] += cross;
                });
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template<class CommType, typename SpinDensityType, typename MagneticField, class occupations_array_type, class kpin_type>
    void eval_psi_vz_psi(CommType & comm, SpinDensityType const & spin_density, MagneticField const & B, occupations_array_type const & occupations, kpin_type const & kpin, double & zeeman_ener) {
        
        basis::field_set<basis::real_space, double> vz(spin_density.skeleton());
        vz.fill(0.0);
        compute_vz(B, vz);

        basis::field<basis::real_space, double> rfield(vz.basis());
        rfield.fill(0.0);
        int iphi = 0;
        for (auto & phi : kpin) {
            compute_psi_vz_psi_ofr(occupations[iphi], phi, vz, rfield);
            iphi++;
        }

        rfield.all_reduce(comm);
        zeeman_ener += operations::integral(rfield);
    }
};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ZEEMAN_COUPLING_UNIT_TEST
#undef INQ_HAMILTONIAN_ZEEMAN_COUPLING_UNIT_TEST

#include <perturbations/magnetic.hpp>
#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

    using namespace inq;
    using namespace inq::magnitude;
    using Catch::Approx;

    parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

    SECTION("Spin polarized zeeman calculation") {
        auto par = input::parallelization(comm);
        auto ions = systems::ions(systems::cell::cubic(10.0_b));
        ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
        auto electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
        ground_state::initial_guess(ions, electrons);
        perturbations::magnetic B{{0.0, 0.0, -1.0}};
        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(200).mixing(0.1), B);
        auto mag = observables::total_magnetization(electrons.spin_density());
        CHECK(mag[0]/mag.length()   == 0.0);
        CHECK(mag[1]/mag.length()   == 0.0);
        CHECK(mag[2]/mag.length()   ==-1.0);
        auto zeeman_ener = result.energy.zeeman_energy();
        Approx target = Approx(zeeman_ener).epsilon(1.e-10);

        hamiltonian::zeeman_coupling zc_(electrons.states().num_density_components());
        basis::field<basis::real_space, vector3<double>> Bfield(electrons.spin_density().basis());
        Bfield.fill(vector3 {0.0, 0.0, 0.0});
        B.magnetic_field(0.0, Bfield);
        auto zeeman_ener2 = 0.0;
        zc_.eval_psi_vz_psi(electrons.kpin_states_comm(), electrons.spin_density(), Bfield, electrons.occupations(), electrons.kpin(), zeeman_ener2);
        CHECK(zeeman_ener2 == target);
    }

    SECTION("Spin non collinear zeeman calculation") {
        auto par = input::parallelization(comm);
        auto ions = systems::ions(systems::cell::cubic(10.0_b));
        ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
        auto electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_non_collinear());
        ground_state::initial_guess(ions, electrons);
        perturbations::magnetic B{{0.0, 0.0, -1.0}};

        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(200).mixing(0.1), B);
        auto mag = observables::total_magnetization(electrons.spin_density());
        auto mx = mag[0]/mag.length();
        auto my = mag[1]/mag.length();
        auto mz = mag[2]/mag.length();
        CHECK(abs(mx) < 1.e-7);
        CHECK(abs(my) < 1.e-7);
        CHECK(abs(mz + 1.0) < 1.e-7);

        auto zeeman_ener = result.energy.zeeman_energy();
        Approx target = Approx(zeeman_ener).epsilon(1.e-10);
        hamiltonian::zeeman_coupling zc_(electrons.states().num_density_components());
        basis::field<basis::real_space, vector3<double>> Bfield(electrons.spin_density().basis());
        Bfield.fill(vector3 {0.0, 0.0, 0.0});
        B.magnetic_field(0.0, Bfield);
        auto zeeman_ener2 = 0.0;
        zc_.eval_psi_vz_psi(electrons.kpin_states_comm(), electrons.spin_density(), Bfield, electrons.occupations(), electrons.kpin(), zeeman_ener2);
        CHECK(zeeman_ener2 == target);
    }
}
#endif