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
        std::cout << "SPIN COMPONENTS : " << spin_components_ << std::endl;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template<typename SpinDensityType, typename MagneticField, typename VKSType>
    void operator()(SpinDensityType const & spin_density, MagneticField const & B, VKSType & vks, double & nvz) const {

        basis::field_set<basis::real_space, double> vz(vks.skeleton());
        vz.fill(0.0);

        assert(vz.set_size() == spin_components_);
    }

};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ZEEMAN_COUPLING_UNIT_TEST
#undef INQ_HAMILTONIAN_ZEEMAN_COUPLING_UNIT_TEST

#include <perturbations/magnetic.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

    using namespace inq;
    using namespace inq::magnitude;

    parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

    SECTION("Spin polarized zeeman calculation") {
        auto par = input::parallelization(comm);
        auto ions = systems::ions(systems::cell::cubic(10.0_b));
        ions.insert("H", {0.0_b, 0.0_b, 0.0_b});
        auto electrons = systems::electrons(par, ions, options::electrons{}.cutoff(30.0_Ha).extra_states(2).spin_polarized());
        ground_state::initial_guess(ions, electrons);
        perturbations::magnetic B{{0.0, 0.0, -1.0}};
        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), B);
    }
}
#endif