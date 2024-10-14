/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__MAGNETIC_FIELD
#define INQ__PERTURBATIONS__MAGNETIC_FIELD

#include <inq_config.h>

namespace inq {
namespace perturbations {

class magnetic_field : public none {

    vector3<double> uniform_magnetic_field_;
    
public:

    magnetic_field(vector3<double> const & field):
        uniform_magnetic_field_(field)
    {
    }

    auto has_uniform_magnetic_field() const {
        return true;
    }

    auto uniform_magnetic_field(double /*time*/) const {
        return uniform_magnetic_field_;
    }

    template<class OStream>
    friend OStream & operator<<(OStream & out, magnetic_field const & self){
        return out;
    }

};

}
}
#endif

#ifdef INQ__PERTURBATIONS__MAGNETIC_FIELD_UNIT_TEST
#undef INQ__PERTURBATIONS__MAGNETIC_FIELD_UNIT_TEST

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

    perturbations::magnetic_field bfield{{0.0, 0.0, 1.0}};

    CHECK(bfield.has_magnetic_field());
    CHECK(bfield.uniform_magnetic_field(/*time*/ 0.0)     == vector3<double>{0.0, 0.0, 1.0});
    CHECK(bfield.uniform_magnetic_field(/*time*/ 1000.0)  == vector3<double>{0.0, 0.0, 1.0});
}
#endif