/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__PARALLELIZATION
#define INQ__INPUT__PARALLELIZATION

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <parallel/communicator.hpp>

#include <parallel/partition.hpp>
#include <utils/factors.hpp>

namespace inq {
namespace input {

  class parallelization {

  public:

		static auto dimension_kpoints(){
			return 2;
		}
		
		static auto dimension_domains(){
			return 0;
		}
		
		static auto dimension_states(){
			return 1;
		}

		static auto optimal_nprocs(int size, int max_comm_size, double threshold){
			for(utils::factors_reverse fac(max_comm_size); fac != fac.end(); ++fac){
				parallel::partition part(size, *fac, 0);

				if(part.local_size(*fac - 1) == 0) continue; //avoid empty partitions
				if(part.waste() <= threshold) return *fac;
			}
			return 1;
		}

		template <class CommType>
    explicit parallelization(CommType & comm):
			nproc_kpts_(boost::mpi3::fill),
			nproc_states_(1),
			nproc_domains_(boost::mpi3::fill),
      comm_(comm)			
		{
    }

		auto & comm() const {
			return comm_;
		}

		auto cart_comm(int nspin, int nkpoints) const {
			assert(nspin == 1 or nspin == 2);
			
			auto nproc_kpts = optimal_nprocs(nkpoints*nspin, comm_.size(), kpoint_efficiency_threshold);
			if(nproc_kpts_ != boost::mpi3::fill) nproc_kpts = nproc_kpts_;

			std::array<int, 3> nprocs;
			nprocs[dimension_kpoints()] = nproc_kpts;
			nprocs[dimension_domains()] = nproc_domains_;
			nprocs[dimension_states()] = nproc_states_;
			
			return parallel::cartesian_communicator<3>(comm_, nprocs);
    }

		auto states(int num = boost::mpi3::fill){
			auto ret = *this;
			ret.nproc_states_ = num;
			return ret;
		}
		
		auto domains(int num = boost::mpi3::fill){
			auto ret = *this;
			ret.nproc_domains_ = num;
			return ret;
		}

		auto kpoints(int num = boost::mpi3::fill){
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
    mutable parallel::communicator comm_;
		constexpr static double const kpoint_efficiency_threshold = 0.1;
  };
    
}
}
#endif

#ifdef INQ_INPUT_PARALLELIZATION_UNIT_TEST
#undef INQ_INPUT_PARALLELIZATION_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
	using namespace Catch::literals;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	input::parallelization par(comm);
	
	auto cart_comm = par.kpoints(1).states(comm.size()).domains(1).cart_comm(2, 10);

	CHECK(par.comm() == comm);
	CHECK(cart_comm.size() == comm.size());

	SECTION("optimize parallelization"){
		CHECK(input::parallelization::optimal_nprocs(16, 4, 0.05) == 4);
		CHECK(input::parallelization::optimal_nprocs(15, 8, 0.1) == 8);
		CHECK(input::parallelization::optimal_nprocs(31, 4, 0.05) == 4);
		CHECK(input::parallelization::optimal_nprocs(20, 38, 0.05) == 2);
		CHECK(input::parallelization::optimal_nprocs(34785, 78, 0.05) == 78);
		CHECK(input::parallelization::optimal_nprocs(12, 64, 0.05) == 4);
		CHECK(input::parallelization::optimal_nprocs(1, 737730, 0.05) == 1);		
	}
	
}
#endif
