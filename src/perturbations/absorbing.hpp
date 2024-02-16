/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__ABSORBING
#define INQ__PERTURBATIONS__ABSORBING

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Yifan Yao
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>
#include <perturbations/none.hpp>

namespace inq {
namespace perturbations {

class absorbing : public perturbations::none {

public:
	absorbing(quantity<magnitude::energy> amplitude, double mid_pos, double width):
		amplitude_(amplitude.in_atomic_units()),
		mid_pos_(mid_pos),
		width_(width)
	{
		mid_pos_ = fmod(mid_pos_, 1); //move mid_pos_ to the first unit cell
		if (mid_pos_ >= 0.5) mid_pos_ = mid_pos_ - 1;
		assert(mid_pos_ >= -0.5 and mid_pos_ < 0.5);
	}

	auto has_potential() const {
		return true;
	}
	
	template<typename PotentialType>
	void potential(const double time, PotentialType & potential) const {
		
		gpu::run(potential.basis().local_sizes()[2], potential.basis().local_sizes()[1], potential.basis().local_sizes()[0],
						 [point_op = potential.basis().point_op(), vk = begin(potential.cubic()), mid_pos = mid_pos_, width = width_, amplitude = amplitude_] GPU_LAMBDA (auto iz, auto iy, auto ix) {
							 auto rr = point_op.rvector(ix, iy, iz);
							 if (rr[2] > mid_pos - width/2 and rr[2] < mid_pos + width/2) {
								 vk[ix][iy][iz] += complex(0.0, amplitude*pow(sin((rr[2] - (mid_pos - width/2))*M_PI/2/(width/2)), 2));
							 }
						 });
	}

private:
	double amplitude_;
	double mid_pos_;
	double width_;
};

}
}
#endif

#ifdef INQ_PERTURBATIONS_ABSORBING_UNIT_TEST
#undef INQ_PERTURBATIONS_ABSORBING_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
    perturbations::absorbing nop(1.0_Ha, 0.3, 0.1);
    CHECK(not nop.has_uniform_electric_field());
}
#endif
