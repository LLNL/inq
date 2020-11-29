/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__RESULT
#define INQ__GROUND_STATE__RESULT

#include <hamiltonian/energy.hpp>
#include <math/array.hpp>
#include <math/vector3.hpp>

namespace inq {
namespace ground_state {

class result {
public:
	hamiltonian::energy energy;
	math::vector3<double> dipole;
	math::array<math::vector3<double>, 1> forces;
};

}
}

#ifdef INQ_GROUND_STATE_RESULT_UNIT_TEST
#undef INQ_GROUND_STATE_RESULT_UNIT_TEST
#endif

#endif

