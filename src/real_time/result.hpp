/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa

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
		std::vector<math::vector3<double>> dipole;
		std::vector<systems::ions> ions;
		std::vector<std::vector<math::vector3<double>>> coordinates;
		std::vector<std::vector<math::vector3<double>>> velocities;		
		std::vector<math::array<math::vector3<double>, 1>> forces;
	};

}
}

#ifdef INQ_REAL_TIME_RESULT_UNIT_TEST
#undef INQ_REAL_TIME_RESULT_UNIT_TEST
#endif

#endif

