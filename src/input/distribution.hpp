/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__DISTRIBUTION
#define INQ__INPUT__DISTRIBUTION

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/communicator.hpp>

namespace inq {
namespace input {

  class distribution {

  public:

    explicit distribution(boost::mpi3::communicator & comm):
			nproc_kpts_(1),
			nproc_states_(1),
			nproc_domains_(boost::mpi3::fill),
      comm_(comm)			
		{
    }

		auto cart_comm() const {
			return boost::mpi3::cartesian_communicator<3>(comm_, {nproc_kpts_, nproc_states_, nproc_domains_});
    }

		auto states(int num){
			auto ret = *this;
			ret.nproc_states_ = num;
			return ret;
		}
		
		auto domains(int num){
			auto ret = *this;
			ret.nproc_domains_ = num;
			return ret;
		}

		auto kpoints(int num){
			auto ret = *this;
			ret.nproc_kpts_ = num;
			return ret;
		}

		auto size() const {
			return comm_.size();
		}

	private:

		int nproc_kpts_;
		int nproc_states_;
		int nproc_domains_;

    mutable boost::mpi3::communicator comm_;
    
  };
    
}
}

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_DISTRIBUTION_UNIT_TEST
#undef INQ_INPUT_DISTRIBUTION_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("class input::distribution", "[inq::input::distribution]") {

  using namespace inq;
	using namespace Catch::literals;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	input::distribution dist(comm);
	
	auto cart_comm = dist.kpoints(1).states(comm.size()).domains(1).cart_comm();

	CHECK(cart_comm.size() == comm.size());

}

#endif
   
#endif
