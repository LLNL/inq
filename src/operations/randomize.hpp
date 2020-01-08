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


#ifdef HAVE_CUDA

		uint64_t start = phi.dist().start();
		uint64_t set_size = phi.set_size();
		uint64_t size_z = phi.basis().sizes()[2];
		uint64_t size_y = phi.basis().sizes()[1];
		auto phicub = begin(phi.cubic());

		gpu::run(phi.dist().local_size(), phi.basis().sizes()[2], phi.basis().sizes()[1], phi.basis().sizes()[0],
						 [=] __device__ (uint64_t ist, uint64_t iz, uint64_t iy, uint64_t ix){
						uniform_distribution<typename field_set_type::element_type> dist;
						uint64_t step = ist + start + set_size*(iz + size_z*(iy + ix*size_y));
						pcg32 rng(seed);
						rng.discard(step*dist.rngs_per_sample);
						phicub[ix][iy][iz][ist] = dist(rng);
					});
#else
		uniform_distribution<typename field_set_type::element_type> dist;

		for(uint64_t ix = 0; ix < uint64_t(phi.basis().sizes()[0]); ix++){
			for(uint64_t iy = 0; iy < uint64_t(phi.basis().sizes()[1]); iy++){
				for(uint64_t iz = 0; iz < uint64_t(phi.basis().sizes()[2]); iz++){
					for(uint64_t ist = 0; ist < uint64_t(phi.dist().local_size()); ist++) {

						uint64_t step = ist + phi.dist().start() + phi.set_size()*(iz + phi.basis().sizes()[2]*(iy + ix*phi.basis().sizes()[1]));
						pcg32 rng(seed);
						rng.discard(step*dist.rngs_per_sample);
						phi.cubic()[ix][iy][iz][ist] = dist(rng);
						
					}
				}
			}
		}
#endif
		
  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::randomize", "[operations::randomize]") {

	using namespace Catch::literals;
  using math::d3vector;

	const int nst = 12;

	double ll = 10.0;
	
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::real_space bas(cell, input::basis::cutoff_energy(20.0));
	
	SECTION("double"){
		basis::field_set<basis::real_space, double> aa(bas, nst, boost::mpi3::environment::get_world_instance());

		aa = 0.0;
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		/*for(int ist = 0; ist < aa.dist().local_size(); ist++){
			std::cout << norms[ist] << std::endl;
			}*/

		if(aa.dist().contains(0))  REQUIRE(norms[0  - aa.dist().start()] == 336.099_a);
		if(aa.dist().contains(1))  REQUIRE(norms[1  - aa.dist().start()] == 335.697_a);
		if(aa.dist().contains(2))  REQUIRE(norms[2  - aa.dist().start()] == 335.101_a);
		if(aa.dist().contains(3))  REQUIRE(norms[3  - aa.dist().start()] == 327.385_a);
		if(aa.dist().contains(4))  REQUIRE(norms[4  - aa.dist().start()] == 337.327_a);
		if(aa.dist().contains(5))  REQUIRE(norms[5  - aa.dist().start()] == 330.692_a);
		if(aa.dist().contains(6))  REQUIRE(norms[6  - aa.dist().start()] == 331.003_a);
		if(aa.dist().contains(7))  REQUIRE(norms[7  - aa.dist().start()] == 328.333_a);
		if(aa.dist().contains(8))  REQUIRE(norms[8  - aa.dist().start()] == 333.662_a);
		if(aa.dist().contains(9))  REQUIRE(norms[9  - aa.dist().start()] == 330.545_a);
		if(aa.dist().contains(10)) REQUIRE(norms[10 - aa.dist().start()] == 335.836_a);
		if(aa.dist().contains(11)) REQUIRE(norms[11 - aa.dist().start()] == 328.899_a);

	}
	
	SECTION("complex"){
		
		basis::field_set<basis::real_space, complex> aa(bas, nst, boost::mpi3::environment::get_world_instance());

		aa = 0.0;
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		/*		for(int ist = 0; ist < aa.dist().local_size(); ist++){
			std::cout << std::scientific << real(norms[ist])<< std::endl;
			}*/

		if(aa.dist().contains(0))  REQUIRE(real(norms[0  - aa.dist().start()]) == 670.4340_a);
		if(aa.dist().contains(1))  REQUIRE(real(norms[1  - aa.dist().start()]) == 663.3693_a);
		if(aa.dist().contains(2))  REQUIRE(real(norms[2  - aa.dist().start()]) == 665.3004_a);
		if(aa.dist().contains(3))  REQUIRE(real(norms[3  - aa.dist().start()]) == 660.0291_a);
		if(aa.dist().contains(4))  REQUIRE(real(norms[4  - aa.dist().start()]) == 660.9823_a);
		if(aa.dist().contains(5))  REQUIRE(real(norms[5  - aa.dist().start()]) == 659.2983_a);
		if(aa.dist().contains(6))  REQUIRE(real(norms[6  - aa.dist().start()]) == 664.7990_a);
		if(aa.dist().contains(7))  REQUIRE(real(norms[7  - aa.dist().start()]) == 666.0472_a);
		if(aa.dist().contains(8))  REQUIRE(real(norms[8  - aa.dist().start()]) == 669.8478_a);
		if(aa.dist().contains(9))  REQUIRE(real(norms[9  - aa.dist().start()]) == 667.2162_a);
		if(aa.dist().contains(10)) REQUIRE(real(norms[10 - aa.dist().start()]) == 666.8721_a);
		if(aa.dist().contains(11)) REQUIRE(real(norms[11 - aa.dist().start()]) == 668.4646_a);

	}

	
}

#endif
#endif



