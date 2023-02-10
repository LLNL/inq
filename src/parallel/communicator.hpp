/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__COMMUNICATOR
#define INQ__PARALLEL__COMMUNICATOR

/*
 Copyright (C) 2022 Xavier Andrade

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

#include <inq_config.h>

#include <mpi3/communicator.hpp>
#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/environment.hpp>

#ifdef ENABLE_NCCL
#include <mpi3/nccl/communicator.hpp>
#endif

#include <cassert>
#include <optional>

namespace inq{
namespace parallel {


template <class CommType>
class hybrid_communicator : public CommType {

#ifdef ENABLE_NCCL
	std::optional<boost::mpi3::nccl::communicator> nccl_comm_;
#endif
	
public:

	using CommType::CommType;
	
	hybrid_communicator(hybrid_communicator const &) = delete;

	hybrid_communicator():
		CommType()
	{
	}

  hybrid_communicator(hybrid_communicator & arg):
    CommType(arg)
  {
  }
	
	template <class ArgType>
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
		if(nccl_comm_.has_value()) return;
		nccl_comm_.emplace(*this);
		assert(nccl_comm_.has_value());
		assert(nccl_comm_->size() == this->size());
#endif
	}

#ifdef ENABLE_NCCL
	auto & nccl_comm() {
		assert(nccl_comm_.has_value());
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

TEST_CASE("class parallel::communicator", "[parallel::communicator]") {
  
}
#endif
