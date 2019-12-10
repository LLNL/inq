/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__RANDOMIZE
#define OPERATIONS__RANDOMIZE

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

namespace operations {

	template <class field_set_type>
  void randomize(field_set_type & phi){

		auto seed = phi.basis().size()*phi.set_size();

		pcg32 rng2(seed);
		double maxrng = std::numeric_limits<decltype(rng2())>::max();

#ifdef HAVE_CUDA

		uint64_t start = phi.dist().start();
		uint64_t set_size = phi.set_size();
		uint64_t size_z = phi.basis().sizes()[2];
		uint64_t size_y = phi.basis().sizes()[1];
		auto phicub = begin(phi.cubic());

		gpu::run(phi.dist().local_size(), phi.basis().sizes()[2], phi.basis().sizes()[1], phi.basis().sizes()[0],
						 [=] __device__ (uint64_t ist, uint64_t iz, uint64_t iy, uint64_t ix){
							 uint64_t step = ist + start + set_size*(iz + size_z*(iy + ix*size_y));
							 pcg32 rng(seed);
							 rng.advance(step);
							 phicub[ix][iy][iz][ist] = rng()/maxrng;
							 });

#else

		for(uint64_t ix = 0; ix < uint64_t(phi.basis().sizes()[0]); ix++){
			for(uint64_t iy = 0; iy < uint64_t(phi.basis().sizes()[1]); iy++){
				for(uint64_t iz = 0; iz < uint64_t(phi.basis().sizes()[2]); iz++){
					for(uint64_t ist = 0; ist < uint64_t(phi.dist().local_size()); ist++) {

						uint64_t step = ist + phi.dist().start() + phi.set_size()*(iz + phi.basis().sizes()[2]*(iy + ix*phi.basis().sizes()[1]));
						pcg32 rng(seed);
						rng.advance(step);
						phi.cubic()[ix][iy][iz][ist] = rng()/maxrng;
						
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
		
		basis::field_set<basis::real_space, double> aa(bas, nst);

		aa = 0.0;
		
		operations::randomize(aa);
		
		auto norms = operations::overlap_diagonal(aa);

		for(int ist = 0; ist < nst; ist++){
			std::cout << norms[ist] << std::endl;
		}

		REQUIRE(norms[0]  == 336.099_a);
		REQUIRE(norms[1]  == 335.697_a);
		REQUIRE(norms[2]  == 335.101_a);
		REQUIRE(norms[3]  == 327.385_a);
		REQUIRE(norms[4]  == 337.327_a);
		REQUIRE(norms[5]  == 330.692_a);
		REQUIRE(norms[6]  == 331.003_a);
		REQUIRE(norms[7]  == 328.333_a);
		REQUIRE(norms[8]  == 333.662_a);
		REQUIRE(norms[9]  == 330.545_a);
		REQUIRE(norms[10] == 335.836_a);
		REQUIRE(norms[11] == 328.899_a);
		
	}
	
	
}

#endif
#endif



