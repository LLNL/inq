/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__RESULT
#define INQ__GROUND_STATE__RESULT

#include <math/vec3d.hpp>

namespace inq {
namespace ground_state {

class result {
public:
	hamiltonian::energy energy;
	math::vec3d dipole;
};

}
}
#endif

