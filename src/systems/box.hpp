/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SYSTEMS__CELL
#define INQ__SYSTEMS__CELL

/*
 Copyright (C) 2019 Xavier Andrade

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

#include <magnitude/length.hpp>
#include <math/vector3.hpp>
#include <utils/merge_optional.hpp>

#include <optional>
#include <cassert>
#include <array>

namespace inq {
namespace systems {

class box {
 
public:

	static auto cubic(quantity<magnitude::length> lat_par){
		auto aa = lat_par.in_atomic_units();
		return box(math::vector3<double>(aa, 0.0, 0.0), math::vector3<double>(0.0, aa, 0.0), math::vector3<double>(0.0, 0.0, aa));
	}

	static box orthorhombic(
		quantity<magnitude::length> aa, 
		quantity<magnitude::length> bb, 
		quantity<magnitude::length> cc
	){
		return {
			math::vector3<double>(aa.in_atomic_units(), 0.0, 0.0), 
			math::vector3<double>(0.0, bb.in_atomic_units(), 0.0), 
			math::vector3<double>(0.0, 0.0, cc.in_atomic_units())
		};
	}

	auto & periodic() {
		periodic_dimensions_ = 3;
		return *this;
	}
		
	auto & finite() {
		periodic_dimensions_ = 0;
		return *this;
	}
		
	auto & operator[](const int ii) const {
		return lattice_vectors_[ii].value();
	}

	auto periodic_dimensions_value() const {
		return periodic_dimensions_.value_or(3);
	}

	friend auto operator|(const box & box1, const box & box2){
		using inq::utils::merge_optional;

		box rbox;
		rbox.lattice_vectors_[0]	= merge_optional(box1.lattice_vectors_[0], box2.lattice_vectors_[0]);
		rbox.lattice_vectors_[1]	= merge_optional(box1.lattice_vectors_[1], box2.lattice_vectors_[1]);
		rbox.lattice_vectors_[2]	= merge_optional(box1.lattice_vectors_[2], box2.lattice_vectors_[2]);
		rbox.periodic_dimensions_	= merge_optional(box1.periodic_dimensions_, box2.periodic_dimensions_);
		return rbox;
	}
	
private:

	box(const math::vector3<double> & a0, const math::vector3<double> & a1, const math::vector3<double> & a2){
		lattice_vectors_[0] = a0;
		lattice_vectors_[1] = a1;
		lattice_vectors_[2] = a2;
	}

	box(){
	}

	std::array<std::optional<math::vector3<double>>, 3> lattice_vectors_;
	std::optional<int> periodic_dimensions_;
		
};

}
}

#ifdef INQ_SYSTEMS_BOX_UNIT_TEST
#undef INQ_SYSTEMS_BOX_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class systems::box", "[systems::box]") {
  
	using namespace inq;
	using namespace magnitude;
	using namespace Catch::literals;

	SECTION("Cubic"){

		auto ci = systems::box::cubic(10.2_b);

		CHECK(ci[0][0] == 10.2_a);
		CHECK(ci[0][1] == 0.0_a);
		CHECK(ci[0][2] == 0.0_a);
		CHECK(ci[1][0] == 0.0_a);
		CHECK(ci[1][1] == 10.2_a);
		CHECK(ci[1][2] == 0.0_a);
		CHECK(ci[2][0] == 0.0_a);
		CHECK(ci[2][1] == 0.0_a);
		CHECK(ci[2][2] == 10.2_a);

		CHECK(ci.periodic_dimensions_value() == 3);
		
	}
	
	SECTION("Cubic finite"){

		auto ci = systems::box::cubic(10.2_b).finite();

		CHECK(ci[0][0] == 10.2_a);
		CHECK(ci[0][1] == 0.0_a);
		CHECK(ci[0][2] == 0.0_a);
		CHECK(ci[1][0] == 0.0_a);
		CHECK(ci[1][1] == 10.2_a);
		CHECK(ci[1][2] == 0.0_a);
		CHECK(ci[2][0] == 0.0_a);
		CHECK(ci[2][1] == 0.0_a);
		CHECK(ci[2][2] == 10.2_a);

		CHECK(ci.periodic_dimensions_value() == 0);
		
	}
	
	SECTION("Parallelepipedic"){

		auto ci = systems::box::orthorhombic(10.2_b, 5.7_b, 8.3_b).periodic();

		CHECK(ci[0][0] == 10.2_a);
		CHECK(ci[0][1] == 0.0_a);
		CHECK(ci[0][2] == 0.0_a);
		CHECK(ci[1][0] == 0.0_a);
		CHECK(ci[1][1] == 5.7_a);
		CHECK(ci[1][2] == 0.0_a);
		CHECK(ci[2][0] == 0.0_a);
		CHECK(ci[2][1] == 0.0_a);
		CHECK(ci[2][2] == 8.3_a);

		CHECK(ci.periodic_dimensions_value() == 3);
		
	}
  
}
#endif

    
#endif
