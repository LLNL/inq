/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__GLOBAL_INDEX
#define INQ__PARALLEL__GLOBAL_INDEX

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

class global_index {

 public:

	constexpr explicit global_index(long val):
		value_(val){
		}

	constexpr auto & value() const {
		return value_;
	}
	
 private:
	long value_; 
};

}
}
#endif

#ifdef INQ_PARALLEL_GLOBAL_INDEX_UNIT_TEST
#undef INQ_PARALLEL_GLOBAL_INDEX_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;

	parallel::global_index gi(10);

	CHECK(gi.value() == 10);
	
}
#endif
