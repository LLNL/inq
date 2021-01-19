/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2020-2021 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__REAL_TIME__RESULT
#define INQ__REAL_TIME__RESULT

#include <math/vector3.hpp>
#include <operations/integral.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>
#include <observables/dipole.hpp>

#include <vector>

namespace inq {
namespace real_time {

	class result {
	public:
		std::vector<double> time;
		std::vector<double> energy;
		std::vector<double> electron_number;		
		std::vector<math::vector3<double>> dipole;
		std::vector<std::vector<math::vector3<double>>> coordinates;
		std::vector<std::vector<math::vector3<double>>> velocities;
		std::vector<math::array<math::vector3<double>, 1>> forces;

		template <class ForcesType>
		void save_iteration_results(double iter_time, systems::ions const & ions, systems::electrons const & electrons, hamiltonian::energy const & iter_energy, ForcesType const & iter_forces){
			time.push_back(iter_time);
			energy.push_back(iter_energy.total());
			electron_number.push_back(operations::integral(electrons.density_));
			dipole.push_back(observables::dipole(ions, electrons));
			coordinates.push_back(ions.geo().coordinates());
			velocities.push_back(ions.geo().velocities());
			forces.push_back(iter_forces);
		}

	};

}
}

#ifdef INQ_REAL_TIME_RESULT_UNIT_TEST
#undef INQ_REAL_TIME_RESULT_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("real_time::result", "[real_time::result]") {

	inq::real_time::result res;
	
}

#endif

#endif

