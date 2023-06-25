/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__RANDOMIZE
#define OPERATIONS__RANDOMIZE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include <pcg-cpp/pcg_random.hpp>

#include <basis/field_set.hpp>
#include <math/complex.hpp>
#include <operations/overlap.hpp>

template<class T>
struct uniform_distribution;

template<>
struct uniform_distribution<double>{
	template<class Generator>
	auto operator()(Generator& g) GPU_FUNCTION {
		static double const max = std::numeric_limits<typename Generator::result_type>::max() + 1.;
		return 2.0*g()/max - 1.0;
	}
	static constexpr std::size_t rngs_per_sample = 1;
};

template<>
struct uniform_distribution<inq::complex>{
	using result_type = inq::complex;
	using param_type = void;
	uniform_distribution<double> impl_;
	template<class Generator> 
	result_type operator()(Generator& g) GPU_FUNCTION {
		return {impl_(g), impl_(g)};
	}
	static constexpr std::size_t rngs_per_sample = 2;
};

namespace inq {
namespace operations {

	template <class field_set_type>
	void randomize(field_set_type & phi, int const index = 0){

		auto seed = phi.basis().size()*phi.set_size() + index;

		uint64_t st_start = phi.set_part().start();
		uint64_t x_start = phi.basis().cubic_part(0).start();
		uint64_t y_start = phi.basis().cubic_part(1).start();
		uint64_t z_start = phi.basis().cubic_part(2).start();
		uint64_t set_size = phi.set_size();
		uint64_t size_z = phi.basis().sizes()[2];
		uint64_t size_y = phi.basis().sizes()[1];
		auto phicub = begin(phi.hypercubic());
		
		gpu::run(phi.set_part().local_size(), phi.basis().cubic_part(2).local_size(), phi.basis().cubic_part(1).local_size(), phi.basis().cubic_part(0).local_size(),
						 [=] GPU_LAMBDA (uint64_t ist, uint64_t iz, uint64_t iy, uint64_t ix){
						uniform_distribution<typename field_set_type::element_type> dist;
						uint64_t step = ist + st_start + set_size*(iz + z_start + size_z*(iy + y_start + (ix + x_start)*size_y));
						pcg32 rng(seed);
						rng.discard(step*dist.rngs_per_sample);
						phicub[ix][iy][iz][ist] = dist(rng);
					});
  }

}
}
#endif

#ifdef INQ_OPERATIONS_RANDOMIZE_UNIT_TEST
#undef INQ_OPERATIONS_RANDOMIZE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <operations/overlap_diagonal.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
  
	const int nst = 12;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	parallel::cartesian_communicator<2> cart_comm(comm, {});

	auto basis_comm = basis::basis_subcomm(cart_comm);

	basis::real_space bas(ions::unit_cell::cubic(10.0_b), /*spacing =*/ 0.49672941, basis_comm);
	
	SECTION("double"){
		basis::field_set<basis::real_space, double> aa(bas, nst, cart_comm);

		aa.fill(0.0);
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		if(aa.set_part().contains(0))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(0))] == 330.1381395023_a);
		if(aa.set_part().contains(1))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(1))] == 330.5444105287_a);
		if(aa.set_part().contains(2))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(2))] == 331.5435469092_a);
		if(aa.set_part().contains(3))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(3))] == 340.0201855451_a);
		if(aa.set_part().contains(4))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(4))] == 332.959828166_a);
		if(aa.set_part().contains(5))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(5))] == 333.1967703871_a);
		if(aa.set_part().contains(6))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(6))] == 328.9429576951_a);
		if(aa.set_part().contains(7))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(7))] == 336.8659386972_a);
		if(aa.set_part().contains(8))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(8))] == 331.6228884743_a);
		if(aa.set_part().contains(9))  CHECK(norms[aa.set_part().global_to_local(parallel::global_index(9))] == 334.5374311686_a);
		if(aa.set_part().contains(10)) CHECK(norms[aa.set_part().global_to_local(parallel::global_index(10))] == 334.6847273315_a);
		if(aa.set_part().contains(11)) CHECK(norms[aa.set_part().global_to_local(parallel::global_index(11))] == 334.1846706334_a);

	}
	
	SECTION("complex"){
		
		basis::field_set<basis::real_space, complex> aa(bas, nst, cart_comm);

		aa.fill(0.0);
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		if(aa.set_part().contains(0))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(0))]) == 669.1459385544_a);
		if(aa.set_part().contains(1))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(1))]) == 664.8544595319_a);
		if(aa.set_part().contains(2))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(2))]) == 669.6590362331_a);
		if(aa.set_part().contains(3))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(3))]) == 668.4278672543_a);
		if(aa.set_part().contains(4))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(4))]) == 668.0965017232_a);
		if(aa.set_part().contains(5))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(5))]) == 665.1942561788_a);
		if(aa.set_part().contains(6))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(6))]) == 665.6688584485_a);
		if(aa.set_part().contains(7))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(7))]) == 672.3987040517_a);
		if(aa.set_part().contains(8))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(8))]) == 671.7988968555_a);
		if(aa.set_part().contains(9))  CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(9))]) == 672.6486737858_a);
		if(aa.set_part().contains(10)) CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(10))]) == 664.3526206348_a);
		if(aa.set_part().contains(11)) CHECK(real(norms[aa.set_part().global_to_local(parallel::global_index(11))]) == 674.5678452953_a);

	}

	
}
#endif



