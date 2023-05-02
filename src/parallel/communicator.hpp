/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__COMMUNICATOR
#define INQ__PARALLEL__COMMUNICATOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <mpi3/communicator.hpp>
#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/environment.hpp>
#include <utils/profiling.hpp>

#ifdef ENABLE_NCCL
#include <mpi3/nccl/communicator.hpp>
#endif

#include <cassert>
#include <memory>

namespace inq{
namespace parallel {


template <class CommType>
class hybrid_communicator : public CommType {

#ifdef ENABLE_NCCL
	std::shared_ptr<boost::mpi3::nccl::communicator> nccl_comm_;
#endif
	
public:

	using CommType::CommType;
	
	hybrid_communicator(hybrid_communicator const &) = delete;
	hybrid_communicator(hybrid_communicator &&) = default;
	
	hybrid_communicator():
		CommType()
	{
	}

  hybrid_communicator(hybrid_communicator & arg):
    CommType(arg)
  {
  }
	
	template <class ArgType, class = std::enable_if_t<not std::is_base_of<hybrid_communicator, ArgType>::value>>
  hybrid_communicator(ArgType && arg):
    CommType(std::forward<ArgType>(arg))
  {
  }

	auto operator=(hybrid_communicator const & comm) = delete;

	auto operator=(hybrid_communicator & comm) {
		CommType::operator=(CommType(comm));
	}

	void nccl_init() {
#ifdef ENABLE_NCCL
		if(nccl_comm_) return;

		CALI_CXX_MARK_FUNCTION;
		nccl_comm_ = std::make_shared<boost::mpi3::nccl::communicator>(*this);
		assert(nccl_comm_);
		assert(nccl_comm_->size() == this->size());
#endif
	}

#ifdef ENABLE_NCCL
	auto & nccl_comm() {
		assert(nccl_comm_);
		assert(nccl_comm_->size() == this->size());
		return *nccl_comm_;
	}
#endif
	
};

using communicator = hybrid_communicator<boost::mpi3::communicator>;

template<boost::mpi3::dimensionality_type D = boost::mpi3::dynamic_extent>
using cartesian_communicator = hybrid_communicator<boost::mpi3::cartesian_communicator<D>>;

}
}
#endif

#ifdef INQ_PARALLEL_COMMUNICATOR_UNIT_TEST
#undef INQ_PARALLEL_COMMUNICATOR_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
}
#endif
