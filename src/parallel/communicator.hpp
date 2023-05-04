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

template<boost::mpi3::dimensionality_type D = boost::mpi3::dynamic_extent> class cartesian_communicator;

class communicator : public boost::mpi3::communicator {

#ifdef ENABLE_NCCL
	std::shared_ptr<boost::mpi3::nccl::communicator> nccl_comm_;
#endif
	
public:

	using base_comm = boost::mpi3::communicator;
	
	communicator(communicator const &) = delete;
	communicator(communicator &&) = default;
	
	communicator():
		base_comm() {
	}

  communicator(communicator & arg):
    base_comm(arg)
#ifdef ENABLE_NCCL
		, nccl_comm_(arg.nccl_comm_)
#endif
  {
  }

  explicit communicator(boost::mpi3::communicator && arg):
    base_comm(std::forward<boost::mpi3::communicator>(arg)) {
  }

  explicit communicator(boost::mpi3::communicator & arg):
    base_comm(arg) {
  }

	template <boost::mpi3::dimensionality_type D>
	communicator(boost::mpi3::cartesian_communicator<D> && arg):
    base_comm(std::forward<boost::mpi3::cartesian_communicator<D>>(arg))
  {
  }

	template <boost::mpi3::dimensionality_type D>
	communicator(boost::mpi3::cartesian_communicator<D> & arg):
    base_comm(arg) {
  }

	template <boost::mpi3::dimensionality_type D>
	communicator(cartesian_communicator<D> && arg):
    base_comm(std::forward<boost::mpi3::cartesian_communicator<D>>(arg))
  {
  }

	template <boost::mpi3::dimensionality_type D>
	communicator(cartesian_communicator<D> & arg):
    base_comm(arg) {
  }
	
	auto operator=(communicator const & comm) = delete;

	auto operator=(communicator & comm) {
		base_comm::operator=(base_comm(comm));
#ifdef ENABLE_NCCL
		nccl_comm_ = comm.nccl_comm_;
#endif
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

template<boost::mpi3::dimensionality_type D>
class cartesian_communicator : public boost::mpi3::cartesian_communicator<D> {

#ifdef ENABLE_NCCL
	std::shared_ptr<boost::mpi3::nccl::communicator> nccl_comm_;
#endif
	
public:

	using base_comm = boost::mpi3::cartesian_communicator<D>;
	
	cartesian_communicator():
		base_comm() {
	}

	cartesian_communicator(cartesian_communicator const &) = delete;

	cartesian_communicator(cartesian_communicator &&) = default;

  cartesian_communicator(cartesian_communicator & arg):
    base_comm(arg)
#ifdef ENABLE_NCCL
		, nccl_comm_(arg.nccl_comm_)
#endif
  {
  }

	template <typename ShapeType>
  cartesian_communicator(boost::mpi3::communicator & comm, ShapeType const & shape):
    base_comm(comm, shape) {
  }
	
  cartesian_communicator(boost::mpi3::communicator & comm, std::array<int, D> shape):
    base_comm(comm, shape) {
  }
	
  explicit cartesian_communicator(boost::mpi3::cartesian_communicator<D> && arg):
    base_comm(std::forward<boost::mpi3::cartesian_communicator<D>>(arg)) {
  }

  explicit cartesian_communicator(boost::mpi3::cartesian_communicator<D> & arg):
    base_comm(arg) {
  }

	auto operator=(cartesian_communicator const & comm) = delete;

	auto operator=(cartesian_communicator & comm) {
		base_comm::operator=(base_comm(comm));
#ifdef ENABLE_NCCL
		nccl_comm_ = comm.nccl_comm_;
#endif
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

	auto axis(int d1){
		return cartesian_communicator<1>(base_comm::axis(d1));
	}
	
	auto plane(int d1, int d2){
		return cartesian_communicator<2>(base_comm::plane(d1, d2));
	}
};

}
}
#endif

#ifdef INQ_PARALLEL_COMMUNICATOR_UNIT_TEST
#undef INQ_PARALLEL_COMMUNICATOR_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
}
#endif
