/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__REAL_TIME__RESULT
#define INQ__REAL_TIME__RESULT

#include <hamiltonian/energy.hpp>
#include <math/vector3.hpp>

#include <systems/ions.hpp>

#include <vector>

namespace inq {
namespace real_time {

	class result {
	public:
		std::vector<double> time;
		std::vector<double> energy;
		std::vector<math::vec3d> dipole;
		std::vector<systems::ions> ions;
	};

}
}

#ifdef INQ_REAL_TIME_RESULT_UNIT_TEST
#undef INQ_REAL_TIME_RESULT_UNIT_TEST
#endif

#endif

