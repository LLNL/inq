/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__VELOCITY_VERLET
#define INQ__SOLVERS__VELOCITY_VERLET

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>

namespace inq {
namespace solvers {
namespace velocity_verlet { 

template <typename AccelType, typename VelType, typename PosType>
void propagate_positions(double dt, AccelType const & accel, VelType & velocities, PosType & positions){

	assert(long(accel.size()) == long(velocities.size()));
	assert(long(accel.size()) == long(positions.size()));	

	for(long ii = 0; ii < long(accel.size()); ii++){
		velocities[ii] += 0.5*dt*accel[ii];		
		positions[ii] += dt*velocities[ii];
	}
	
}

template <typename AccelType, typename VelType>
void propagate_velocities(double dt, AccelType const & accel, VelType & velocities){

	assert(long(accel.size()) == long(velocities.size()));

	for(long ii = 0; ii < long(accel.size()); ii++) velocities[ii] += 0.5*dt*accel[ii];
}

}
}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_VELOCITY_VERLET_UNIT_TEST
#undef INQ_SOLVERS_VELOCITY_VERLET_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <gpu/array.hpp>

using namespace inq;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	SECTION("Cosine"){

		double dt = 0.001;
		
		std::vector<double> pos(2);
		std::vector<double> vel(2);		
		std::vector<double> acc(2);		

		pos[0] = 1.0;
		vel[0] = 0.0;
		pos[1] = 0.0;
		vel[1] = 1.0;
		
		acc[0] = -pos[0];
		acc[1] = -pos[1];		

		for(int it = 0; it < 10000; it++){

			solvers::velocity_verlet::propagate_positions(dt, acc, vel, pos);

			acc[0] = -pos[0];
			acc[1] = -pos[1];
			
			solvers::velocity_verlet::propagate_velocities(dt, acc, vel);

			CHECK(fabs(pos[0] - cos((it + 1)*dt)) < 1.0e-5);
			CHECK(fabs(vel[0] + sin((it + 1)*dt)) < 1.0e-5);			

			CHECK(fabs(pos[1] - sin((it + 1)*dt)) < 1.0e-5);
			CHECK(fabs(vel[1] - cos((it + 1)*dt)) < 1.0e-5);	
		}
		
  }
	
}
#endif
