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
	
  auto num() const {
    return product(dims_);
  }

	vector3<int> is_shifted() const {
		if(shifted_) return {1, 1, 1};
		return {0, 0, 0};
	}

};

auto gamma(){
	return grid({1, 1, 1}, false);
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

    CHECK(kpts.dims()[0] == 1);
    CHECK(kpts.dims()[1] == 1);
    CHECK(kpts.dims()[2] == 1);
    CHECK(kpts.is_shifted()[0] == 0);
    CHECK(kpts.is_shifted()[1] == 0);
    CHECK(kpts.is_shifted()[2] == 0);
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
}
#endif
