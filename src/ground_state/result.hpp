/* -*- indent-tabs-mode: t -*- */

#ifndef GROUND_STATE__RESULT
#define GROUND_STATE__RESULT

#include <math/vec3d.hpp>

namespace ground_state {

	class result {
	public:
		hamiltonian::energy energy;
		math::vec3d dipole;
	};
}

#endif

