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

	gpu::array<double, 1> buffer(comm.size());
	MPI_Allgather(&nrm2, 1, MPI_DOUBLE, buffer.data(), 1, MPI_DOUBLE, comm.get());
	nrm2 = boost::multi::blas::nrm2(buffer);

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
  
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	auto size = 1000;
	gpu::array<double, 1> buffer(size);
	
	gpu::run(size, [buf = begin(buffer)] GPU_LAMBDA (auto ii) {
		buf[ii] = pow(cos(ii), 3);
	});

	auto serial_nrm2 = boost::multi::blas::nrm2(buffer);

	CHECK(serial_nrm2 == 17.6975_a);

	auto part = parallel::partition(size, comm);

	auto par_nrm2 = +boost::multi::blas::nrm2(buffer({part.start(), part.end()}));
	parallel::reduce_nrm2(par_nrm2, comm);

	CHECK(par_nrm2 == 17.6974674479_a);
}

#endif
