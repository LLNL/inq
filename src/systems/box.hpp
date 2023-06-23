/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SYSTEMS__CELL
#define INQ__SYSTEMS__CELL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <ions/unit_cell.hpp>
#include <magnitude/length.hpp>
#include <magnitude/energy.hpp>
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
		return box(vector3<double>(aa, 0.0, 0.0), vector3<double>(0.0, aa, 0.0), vector3<double>(0.0, 0.0, aa));
	}

	static box orthorhombic(
		quantity<magnitude::length> aa, 
		quantity<magnitude::length> bb, 
		quantity<magnitude::length> cc
	){
		return {
			vector3<double>(aa.in_atomic_units(), 0.0, 0.0), 
			vector3<double>(0.0, bb.in_atomic_units(), 0.0), 
			vector3<double>(0.0, 0.0, cc.in_atomic_units())
		};
	}
	
	static box lattice(vector3<quantity<magnitude::length>> aa, vector3<quantity<magnitude::length>> bb, vector3<quantity<magnitude::length>> cc){
		return {
			vector3<double>(aa[0].in_atomic_units(), aa[1].in_atomic_units(), aa[2].in_atomic_units()), 
			vector3<double>(bb[0].in_atomic_units(), bb[1].in_atomic_units(), bb[2].in_atomic_units()), 
			vector3<double>(cc[0].in_atomic_units(), cc[1].in_atomic_units(), cc[2].in_atomic_units())
		};
	}

	template <typename LatticeType>
	box(LatticeType const & lat):
		box{lat[0], lat[1], lat[2]}
	{
	}
	
	auto & periodic() {
		periodicity_ = 3;
		return *this;
	}
		
	auto & finite() {
		periodicity_ = 0;
		return *this;
	}

	auto & periodicity(int const pval) {
		if(pval > 3 or pval < 0) throw std::runtime_error("inq error: the requested periodicity (" + std::to_string(pval) + ") does not make sense.");
		if(pval == 1) throw std::runtime_error("inq error: periodicity 1 is not implemented yet.");
		periodicity_ = pval;
		return *this;
	}
	
	auto periodicity_value() const {
		return periodicity_.value_or(3);
	}

	auto spherical_grid(bool arg_sph_grid){
		spherical_grid_ = arg_sph_grid;
		return *this;
	}
	
	auto spherical_grid_value() const {
		return spherical_grid_.value_or(false);
	}

	auto & density_factor(double arg_factor){
		density_factor_ = arg_factor;
		return *this;
	}

	auto density_factor_value() const {
		return density_factor_.value_or(1.0);
	}
	
	friend auto operator==(box const& self, box const& other) {
		return
			    self.lattice_vectors_     == other.lattice_vectors_
			and self.periodicity_ == other.periodicity_
			and self.spherical_grid_      == other.spherical_grid_
			and self.density_factor_      == other.density_factor_
		;
	}
	
	friend auto operator!=(box const& self, box const& other) {return not(self == other);}

	auto cell() const {
		return ions::unit_cell(*lattice_vectors_[0], *lattice_vectors_[1], *lattice_vectors_[2], periodicity_value());
	}
	
private:

	box(const vector3<double> & a0, const vector3<double> & a1, const vector3<double> & a2){
		lattice_vectors_[0] = a0;
		lattice_vectors_[1] = a1;
		lattice_vectors_[2] = a2;
	}

	box(){
	}

	std::array<std::optional<vector3<double>>, 3> lattice_vectors_;
	std::optional<int> periodicity_;
	std::optional<bool> spherical_grid_;
	std::optional<double> density_factor_;
		
};

}
}
#endif

#ifdef INQ_SYSTEMS_BOX_UNIT_TEST
#undef INQ_SYSTEMS_BOX_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace magnitude;
	using namespace Catch::literals;

	SECTION("Cubic"){

		auto ci = systems::box::cubic(10.2_b);

		CHECK(ci.cell()[0][0] == 10.2_a);
		CHECK(ci.cell()[0][1] == 0.0_a);
		CHECK(ci.cell()[0][2] == 0.0_a);
		CHECK(ci.cell()[1][0] == 0.0_a);
		CHECK(ci.cell()[1][1] == 10.2_a);
		CHECK(ci.cell()[1][2] == 0.0_a);
		CHECK(ci.cell()[2][0] == 0.0_a);
		CHECK(ci.cell()[2][1] == 0.0_a);
		CHECK(ci.cell()[2][2] == 10.2_a);

		CHECK(ci.periodicity_value() == 3);
		
	}
	
	SECTION("Cubic finite"){

		auto ci = systems::box::cubic(10.2_b).finite();

		CHECK(ci.cell()[0][0] == 10.2_a);
		CHECK(ci.cell()[0][1] == 0.0_a);
		CHECK(ci.cell()[0][2] == 0.0_a);
		CHECK(ci.cell()[1][0] == 0.0_a);
		CHECK(ci.cell()[1][1] == 10.2_a);
		CHECK(ci.cell()[1][2] == 0.0_a);
		CHECK(ci.cell()[2][0] == 0.0_a);
		CHECK(ci.cell()[2][1] == 0.0_a);
		CHECK(ci.cell()[2][2] == 10.2_a);
		CHECK(ci.periodicity_value() == 0);
		CHECK(not ci.spherical_grid_value());
		
	}
	
	SECTION("Parallelepipedic"){

		auto ci = systems::box::orthorhombic(10.2_b, 5.7_b, 8.3_b).periodic();

		CHECK(ci.cell()[0][0] == 10.2_a);
		CHECK(ci.cell()[0][1] == 0.0_a);
		CHECK(ci.cell()[0][2] == 0.0_a);
		CHECK(ci.cell()[1][0] == 0.0_a);
		CHECK(ci.cell()[1][1] == 5.7_a);
		CHECK(ci.cell()[1][2] == 0.0_a);
		CHECK(ci.cell()[2][0] == 0.0_a);
		CHECK(ci.cell()[2][1] == 0.0_a);
		CHECK(ci.cell()[2][2] == 8.3_a);
		CHECK(ci.periodicity_value() == 3);

	}
			
	SECTION("Spherical grid"){

		auto ci = systems::box::orthorhombic(10.2_b, 5.7_b, 8.3_b).periodic().spherical_grid(true);

		CHECK(ci.cell()[0][0] == 10.2_a);
		CHECK(ci.cell()[0][1] == 0.0_a);
		CHECK(ci.cell()[0][2] == 0.0_a);
		CHECK(ci.cell()[1][0] == 0.0_a);
		CHECK(ci.cell()[1][1] == 5.7_a);
		CHECK(ci.cell()[1][2] == 0.0_a);
		CHECK(ci.cell()[2][0] == 0.0_a);
		CHECK(ci.cell()[2][1] == 0.0_a);
		CHECK(ci.cell()[2][2] == 8.3_a);
		CHECK(ci.periodicity_value() == 3);
		CHECK(ci.spherical_grid_value());

	}
	
	SECTION("Non-orthogonal"){

		auto ci = systems::box::lattice({0.0_A, 1.0_A, 1.0_A}, {1.0_A, 0.0_b, 1.0_A}, {1.0_A, 1.0_A, 0.0_A});

		CHECK(ci.cell()[0][0] == 0.0_a);
		CHECK(ci.cell()[0][1] == 1.8897261246_a);
		CHECK(ci.cell()[0][2] == 1.8897261246_a);
		CHECK(ci.cell()[1][0] == 1.8897261246_a);
		CHECK(ci.cell()[1][1] == 0.0_a);
		CHECK(ci.cell()[1][2] == 1.8897261246_a);
		CHECK(ci.cell()[2][0] == 1.8897261246_a);
		CHECK(ci.cell()[2][1] == 1.8897261246_a);
		CHECK(ci.cell()[2][2] == 0.0_a);
		CHECK(ci.periodicity_value() == 3);

	}
	
	SECTION("Equality"){

		auto ci1 = systems::box::orthorhombic(10.2_b, 5.7_b, 8.3_b).periodic().spherical_grid(true);
		auto ci2 = systems::box::orthorhombic(10.2_b, 5.7_b, 8.3_b).periodic().spherical_grid(true);

		CHECK(ci1 == ci2);
		CHECK( not (ci1 != ci2) );
	}
}
#endif
