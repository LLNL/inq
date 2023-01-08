/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__PARALLELIZATION
#define INQ__INPUT__PARALLELIZATION

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

#include <parallel/communicator.hpp>

#include <parallel/partition.hpp>
#include <utils/factors.hpp>

namespace inq {
namespace input {

  class parallelization {

  public:

		static auto dimension_kpoints(){
			return 0;
		}
		
		static auto dimension_domains(){
			return 1;
		}
		
		static auto dimension_states(){
			return 2;
		}

		static auto optimal_nprocs(int size, int max_comm_size, double threshold){
			for(utils::factors_reverse fac(max_comm_size); fac != fac.end(); ++fac){
				parallel::partition part(size, *fac);

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

		auto cart_comm(int nspin, int nkpoints) const {
			auto nproc_kpts = optimal_nprocs(nkpoints*nspin, comm_.size(), kpoint_efficiency_threshold);
			if(nproc_kpts_ != boost::mpi3::fill) nproc_kpts = nproc_kpts_;

			return parallel::cartesian_communicator<3>(comm_, {nproc_kpts, nproc_domains_, nproc_states_});
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

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_PARALLELIZATION_UNIT_TEST
#undef INQ_INPUT_PARALLELIZATION_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("class input::parallelization", "[inq::input::parallelization]") {

  using namespace inq;
	using namespace Catch::literals;

	parallel::communicator comm = boost::mpi3::environment::get_world_instance();
	
	input::parallelization par(comm);
	
	auto cart_comm = par.kpoints(1).states(comm.size()).domains(1).cart_comm(2, 10);

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
   
#endif
