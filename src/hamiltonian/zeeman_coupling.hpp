/* -*- indent-tabs-mode: t -*- */

#ifndef ZEEMAN_COUPL_HPP
#define ZEEMAN_COUPL_HPP

namespace inq {
namespace hamiltonian {

class zeeman_coupling {

private:

    vector3<double> B_;
    int spin_components_;

public:

    zeeman_coupling(int const spin_components, vector3<double> const & B):
        spin_components_(spin_components),
        B_(B)
    {
        assert(spin_components_ > 1);
        std::cout << "SPIN COMPONENTS : " << spin_components_ << std::endl;
        std::cout << B_ << std::endl;   
    }

};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ZEEMAN_COUPL_UNIT_TEST
#undef INQ_HAMILTONIAN_ZEEMAN_COUPL_UNIT_TEST

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

    parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

    SECTION("Spin polarized vz calculation") {
        auto par = input::parallelization(comm);
    }
}
#endif