/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__REDUCE_NRM2
#define INQ__PARALLEL__REDUCE_NRM2

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <cstdlib>

#include <gpu/array.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <inq_config.h>
#include <mpi.h>
#include <parallel/communicator.hpp>
#include <gpu/run.hpp>

namespace inq {
namespace parallel {

template <class Comm>
void reduce_nrm2(double & nrm2, Comm & comm){
	CALI_CXX_MARK_FUNCTION;

	if(comm.size() == 1) return;
	
	auto square = nrm2*nrm2;
	comm.all_reduce_in_place_n(&square, 1, std::plus<>{});
	nrm2 = sqrt(square);
}

}
}
#endif

#ifdef INQ_PARALLEL_REDUCE_NRM2_UNIT_TEST
#undef INQ_PARALLEL_REDUCE_NRM2_UNIT_TEST

#include <gpu/run.hpp>
#include <mpi3/environment.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

  using namespace inq;
  using namespace Catch::literals;
  
	///	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
  
	//  gpu::array<int, 2> buffer({comm.size(), blocksize}, comm.rank());

}

#endif
