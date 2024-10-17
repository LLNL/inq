/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__MAGNET
#define INQ__PERTURBATIONS__MAGNET

#include <inq_config.h>

namespace inq {
namespace perturbations {

class magnetic : public none {

    vector3<double> magnetic_vector_;
    
public:

    magnetic(vector3<double> const & B):
        magnetic_vector_(B)
    {
    }

    auto has_magnetic_field() const {
        return true;
    }

    template<typename MagneticField>
    void magnetic_field(const double time, MagneticField & magnetic) const {
        gpu::run(magnetic.basis().local_size(),
            [b = begin(magnetic.linear()), mv = magnetic_vector_] GPU_LAMBDA (auto ip){
                b[ip] += mv;
            });
    }

    template<class OStream>
    friend OStream & operator<<(OStream & out, magnetic const & self){
        return out;
    }

};

}
}
#endif

#ifdef INQ_PERTURBATIONS_MAGNETIC_UNIT_TEST
#undef INQ_PERTURBATIONS_MAGNETIC_UNIT_TEST

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

    parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

    perturbations::magnetic B{{0.0, 0.0, 1.0}};

    basis::real_space bas(systems::cell::cubic(5.0_b), /*spacing*/ 0.1, comm);
    basis::field<basis::real_space, vector3<double>> Bfield(bas);
    Bfield.fill(vector3<double> {0.0, 0.0, 0.0});

    CHECK(B.has_magnetic_field());
    B.magnetic_field(/*time*/ 0.0, Bfield);
    CHECK(Bfield.linear()[0]     == vector3<double>{0.0, 0.0, 1.0});
    CHECK(Bfield.linear()[1]     == vector3<double>{0.0, 0.0, 1.0});

    B.magnetic_field(/*time*/ 1000.0, Bfield);
    CHECK(Bfield.linear()[0]     == vector3<double>{0.0, 0.0, 2.0});
}
#endif