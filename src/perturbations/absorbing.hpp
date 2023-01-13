/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__ABSORBING
#define INQ__PERTURBATIONS__ABSORBING

/*
 Copyright (C) 2019-2023 Xavier Andrade, Alfredo Correa, Yifan Yao.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>

namespace inq {
namespace perturbations {

class absorbing {

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

    template <typename DummyType>
    void zero_step(DummyType &) const {
    }

    auto has_uniform_electric_field() const {
        return false;
    }

    auto uniform_electric_field(double /*time*/) const {
        return math::vector3<double, math::cartesian> {0.0, 0.0, 0.0};
    }

    auto has_uniform_vector_potential() const {
        return false;
    }

    auto uniform_vector_potential(double /*time*/) const {
        return math::vector3<double, math::cartesian> {0.0, 0.0, 0.0};
    }

    auto has_potential() const {
		return true;
    }

    template<typename PotentialType>
    void potential(const double time, PotentialType & potential) const {
		auto Vcap = [mid_pos = mid_pos_, width = width_, amplitude = amplitude_](inq::math::vector3<double, math::contravariant> rr) {
			if (rr[2] > mid_pos - width/2 && rr[2] < mid_pos + width/2) {
				return complex(0.0, amplitude*pow(sin((rr[2] - (mid_pos - width/2))*M_PI/2/(width/2)),2));
			}
			else
				return complex{0,0};
		};
       	gpu::run(potential.basis().local_sizes()[2], potential.basis().local_sizes()[1], potential.basis().local_sizes()[0],
       	[point_op = potential.basis().point_op(), vk = begin(potential.cubic()), Vcap] GPU_LAMBDA (auto iz, auto iy, auto ix) {
       	    auto rr = point_op.rvector(ix, iy, iz);
       	    vk[ix][iy][iz] += Vcap(rr);
       	});
    }

private:
	double amplitude_;
	double mid_pos_;
	double width_;
};

}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_ABSORBING_UNIT_TEST
#undef INQ_PERTURBATIONS_ABSORBING_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE("perturbations::absorbing", "[perturbations::absorbing]") {
    perturbations::absorbing nop(1.0_Ha, 0.3, 0.1);
    CHECK(not nop.has_uniform_electric_field());
}

#endif
#endif
