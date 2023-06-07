/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__KPOINTS
#define INQ__INPUT__KPOINTS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/vector3.hpp>

#include <vector>
#include <cmath>

namespace inq {
namespace input {
namespace kpoints {

class grid {
	
	vector3<int> dims_;
	bool shifted_;
	
public:
	
  grid(vector3<int> const & dims, bool shifted = false):
    dims_(dims),
		shifted_(shifted){
  }
	
  auto & dims() const {
    return dims_;
  }
	
  auto size() const {
    return product(dims_);
  }

	vector3<int> is_shifted() const {
		if(shifted_) return {1, 1, 1};
		return {0, 0, 0};
	}

};

class list {
	
	std::vector<vector3<double, covariant>> kpoints_;
	std::vector<double> weights_;

public:

	void insert(vector3<double, covariant> const & kpoint, double const & weight = 1.0){
		kpoints_.push_back(kpoint);
		weights_.push_back(weight);
	}

	auto size() const {
		return (long) kpoints_.size();
	}

	auto kpoint(int ik) const {
		return kpoints_[ik];
  }
  
  auto weight(int ik) const {
    return weights_[ik];
  }
	
};

auto single(vector3<double, covariant> const & kpoint, double const & weight = 1.0){
	auto kpts = input::kpoints::list();
	kpts.insert(kpoint, weight);
	return kpts;
}

auto gamma(){
	return single({0.0, 0.0, 0.0}, 1.0);
}


}
}
}
#endif

#ifdef INQ_INPUT_KPOINTS_UNIT_TEST
#undef INQ_INPUT_KPOINTS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
	using namespace Catch::literals;

	SECTION("Gamma - no arguments"){
		auto kpts = input::kpoints::gamma();

		CHECK(kpts.size() == 1);

		CHECK(kpts.kpoint(0)[0] == 0.0_a);
		CHECK(kpts.kpoint(0)[1] == 0.0_a);
		CHECK(kpts.kpoint(0)[2] == 0.0_a);

		CHECK(kpts.weight(0) ==  1.0_a);
	}
  
	SECTION("Grid - one argument"){
		auto kpts = input::kpoints::grid({10, 9, 8});

    CHECK(kpts.dims()[0] == 10);
    CHECK(kpts.dims()[1] == 9);
    CHECK(kpts.dims()[2] == 8);
    CHECK(kpts.is_shifted()[0] == 0);
    CHECK(kpts.is_shifted()[1] == 0);
    CHECK(kpts.is_shifted()[2] == 0);
	}
	
	SECTION("Grid - two arguments"){
		auto kpts = input::kpoints::grid({10, 9, 8}, true);

    CHECK(kpts.dims()[0] == 10);
    CHECK(kpts.dims()[1] == 9);
    CHECK(kpts.dims()[2] == 8);
    CHECK(kpts.is_shifted()[0] == 1);
    CHECK(kpts.is_shifted()[1] == 1);
    CHECK(kpts.is_shifted()[2] == 1);
	}

	SECTION("List"){
		auto kpts = input::kpoints::list();

		kpts.insert({0.5, 0.5, -0.5}, 0.3);
		kpts.insert({0.1, 0.2,  0.3}, 0.7);

		CHECK(kpts.size() == 2);

		CHECK(kpts.kpoint(0)[0] ==  0.5_a);
		CHECK(kpts.kpoint(0)[1] ==  0.5_a);
		CHECK(kpts.kpoint(0)[2] == -0.5_a);
		CHECK(kpts.kpoint(1)[0] ==  0.1_a);
		CHECK(kpts.kpoint(1)[1] ==  0.2_a);
		CHECK(kpts.kpoint(1)[2] ==  0.3_a);

		CHECK(kpts.weight(0) ==  0.3_a);
		CHECK(kpts.weight(1) ==  0.7_a);
		
	}

	
}
#endif
