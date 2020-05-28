/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__REAL_TIME__RESULT
#define INQ__REAL_TIME__RESULT

#include <math/vec3d.hpp>

namespace inq {
namespace real_time {

	class result {
	public:
		std::vector<double> time;
		std::vector<double> energy;
		std::vector<math::vec3d> dipole;
	};

}
}

#endif

