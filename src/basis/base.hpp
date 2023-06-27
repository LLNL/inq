/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__BASE
#define INQ__BASIS__BASE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>

#include <mpi3/environment.hpp>
#include <parallel/communicator.hpp>
#include <parallel/partition.hpp>

namespace inq {
namespace basis {

class base {
	
public:
	
	base(const long size, parallel::communicator & comm):
		comm_(comm),
		part_(size, comm){
	}
	
	auto & part() {
		return part_;
	}
	
	auto & part() const {
		return part_;
	}

	auto & comm() const {
		return comm_;
	}

	auto local_size() const {
		return part_.local_size();
	}

	protected:
	
	mutable parallel::communicator comm_;
	inq::parallel::partition part_;
		
};

}
}
#endif

#ifdef INQ_BASIS_BASE_UNIT_TEST
#undef INQ_BASIS_BASE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;

}
#endif
