/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__RANDOMIZE
#define OPERATIONS__RANDOMIZE

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa.

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

#include <basis/field_set.hpp>
#include <cstdlib>

#include <pcg-cpp/pcg_random.hpp>

template<class T>
struct uniform_distribution;

template<>
struct uniform_distribution<double>{// : std::uniform_real_distribution<double>{
	template<class Generator>
	auto operator()(Generator& g) GPU_FUNCTION {
		static double const max = std::numeric_limits<typename Generator::result_type>::max() + 1.;
		return g()/max;
	}
	static constexpr std::size_t rngs_per_sample = 1;
};

template<>
struct uniform_distribution<std::complex<double>>{
	using result_type = std::complex<double>;
	using param_type = void;
	uniform_distribution<double> impl_;
	template<class Generator> 
	result_type operator()(Generator& g) GPU_FUNCTION {
		return {impl_(g), impl_(g)};
	}
	static constexpr std::size_t rngs_per_sample = 2;
};

namespace operations {

	template <class field_set_type>
	void randomize(field_set_type & phi){

		auto seed = phi.basis().size()*phi.set_size();

		uint64_t st_start = phi.set_dist().start();
		uint64_t x_start = phi.basis().cubic_dist(0).start();
		uint64_t y_start = phi.basis().cubic_dist(1).start();
		uint64_t z_start = phi.basis().cubic_dist(2).start();
		uint64_t set_size = phi.set_size();
		uint64_t size_z = phi.basis().sizes()[2];
		uint64_t size_y = phi.basis().sizes()[1];
		auto phicub = begin(phi.cubic());


		// DATAOPERATIONS GPU:RUN 4D
		gpu::run(phi.set_dist().local_size(), phi.basis().cubic_dist(2).local_size(), phi.basis().cubic_dist(1).local_size(), phi.basis().cubic_dist(0).local_size(),
						 [=] GPU_LAMBDA (uint64_t ist, uint64_t iz, uint64_t iy, uint64_t ix){
						uniform_distribution<typename field_set_type::element_type> dist;
						uint64_t step = ist + st_start + set_size*(iz + z_start + size_z*(iy + y_start + (ix + x_start)*size_y));
						pcg32 rng(seed);
						rng.discard(step*dist.rngs_per_sample);
						phicub[ix][iy][iz][ist] = dist(rng);
					});
	
  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::randomize", "[operations::randomize]") {

	using namespace Catch::literals;
  using math::vec3d;

	const int nst = 12;

	double ll = 10.0;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	boost::mpi3::cartesian_communicator cart_comm(comm, boost::mpi3::dims_create(comm.size(), 2), true);

	auto basis_comm = cart_comm.sub({1, 0});
	
  ions::UnitCell cell(vec3d(ll, 0.0, 0.0), vec3d(0.0, ll, 0.0), vec3d(0.0, 0.0, ll));
  basis::real_space bas(cell, input::basis::cutoff_energy(20.0), basis_comm);
	
	SECTION("double"){
		basis::field_set<basis::real_space, double> aa(bas, nst, cart_comm);

		aa = 0.0;
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		/*for(int ist = 0; ist < aa.set_dist().local_size(); ist++){
			std::cout << norms[ist] << std::endl;
			}*/

		if(aa.set_dist().contains(0))  REQUIRE(norms[aa.set_dist().global_to_local(0)] == 336.099_a);
		if(aa.set_dist().contains(1))  REQUIRE(norms[aa.set_dist().global_to_local(1)] == 335.697_a);
		if(aa.set_dist().contains(2))  REQUIRE(norms[aa.set_dist().global_to_local(2)] == 335.101_a);
		if(aa.set_dist().contains(3))  REQUIRE(norms[aa.set_dist().global_to_local(3)] == 327.385_a);
		if(aa.set_dist().contains(4))  REQUIRE(norms[aa.set_dist().global_to_local(4)] == 337.327_a);
		if(aa.set_dist().contains(5))  REQUIRE(norms[aa.set_dist().global_to_local(5)] == 330.692_a);
		if(aa.set_dist().contains(6))  REQUIRE(norms[aa.set_dist().global_to_local(6)] == 331.003_a);
		if(aa.set_dist().contains(7))  REQUIRE(norms[aa.set_dist().global_to_local(7)] == 328.333_a);
		if(aa.set_dist().contains(8))  REQUIRE(norms[aa.set_dist().global_to_local(8)] == 333.662_a);
		if(aa.set_dist().contains(9))  REQUIRE(norms[aa.set_dist().global_to_local(9)] == 330.545_a);
		if(aa.set_dist().contains(10)) REQUIRE(norms[aa.set_dist().global_to_local(10)] == 335.836_a);
		if(aa.set_dist().contains(11)) REQUIRE(norms[aa.set_dist().global_to_local(11)] == 328.899_a);

	}
	
	SECTION("complex"){
		
		basis::field_set<basis::real_space, complex> aa(bas, nst, cart_comm);

		aa = 0.0;
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		/*		for(int ist = 0; ist < aa.set_dist().local_size(); ist++){
			std::cout << std::scientific << real(norms[ist])<< std::endl;
			}*/

		if(aa.set_dist().contains(0))  REQUIRE(real(norms[aa.set_dist().global_to_local(0)]) == 670.4340_a);
		if(aa.set_dist().contains(1))  REQUIRE(real(norms[aa.set_dist().global_to_local(1)]) == 663.3693_a);
		if(aa.set_dist().contains(2))  REQUIRE(real(norms[aa.set_dist().global_to_local(2)]) == 665.3004_a);
		if(aa.set_dist().contains(3))  REQUIRE(real(norms[aa.set_dist().global_to_local(3)]) == 660.0291_a);
		if(aa.set_dist().contains(4))  REQUIRE(real(norms[aa.set_dist().global_to_local(4)]) == 660.9823_a);
		if(aa.set_dist().contains(5))  REQUIRE(real(norms[aa.set_dist().global_to_local(5)]) == 659.2983_a);
		if(aa.set_dist().contains(6))  REQUIRE(real(norms[aa.set_dist().global_to_local(6)]) == 664.7990_a);
		if(aa.set_dist().contains(7))  REQUIRE(real(norms[aa.set_dist().global_to_local(7)]) == 666.0472_a);
		if(aa.set_dist().contains(8))  REQUIRE(real(norms[aa.set_dist().global_to_local(8)]) == 669.8478_a);
		if(aa.set_dist().contains(9))  REQUIRE(real(norms[aa.set_dist().global_to_local(9)]) == 667.2162_a);
		if(aa.set_dist().contains(10)) REQUIRE(real(norms[aa.set_dist().global_to_local(10)]) == 666.8721_a);
		if(aa.set_dist().contains(11)) REQUIRE(real(norms[aa.set_dist().global_to_local(11)]) == 668.4646_a);

	}

	
}

#endif
#endif



