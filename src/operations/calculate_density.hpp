/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__CALCULATE_DENSITY
#define OPERATIONS__CALCULATE_DENSITY

/*
 Copyright (C) 2019 Xavier Andrade

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

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <math/complex.hpp>
#include <cstdlib>

namespace operations {

  template<class occupations_array_type, class field_set_type>
  basis::field<typename field_set_type::basis_type, double> calculate_density(const occupations_array_type & occupations, field_set_type & phi){

    basis::field<typename field_set_type::basis_type, double> density(phi.basis());

    //DATAOPERATIONS LOOP + GPU::RUN 2D
#ifdef HAVE_CUDA

		const auto nst = phi.dist().local_size();
		auto occupationsp = begin(occupations);
		auto phip = begin(phi.matrix());
		auto densityp = begin(density.linear());
		
		gpu::run(phi.basis().size(),
						 [=] __device__ (auto ipoint){
							 densityp[ipoint] = 0.0;
							 for(int ist = 0; ist < nst; ist++) densityp[ipoint] += occupationsp[ist]*norm(phip[ipoint][ist]);
						 });
		
#else
		
    for(int ipoint = 0; ipoint < phi.basis().size(); ipoint++){
			density.linear()[ipoint] = 0.0;
      for(int ist = 0; ist < phi.dist().local_size(); ist++) density.linear()[ipoint] += occupations[ist]*norm(phi.matrix()[ipoint][ist]);
    }
		
#endif

		if(phi.dist().parallel()){
			phi.dist().comm().all_reduce_in_place_n(static_cast<double *>(density.linear().data()), density.linear().size(), std::plus<>{});
		}
		
    return density;
  }
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::calculate_density", "[operations::calculate_density]") {

	using namespace Catch::literals;

	const int N = 100;
	const int M = 12;
	
	basis::trivial bas(N);
	
	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, M, boost::mpi3::environment::get_world_instance());

		math::array<double, 1> occ(aa.dist().local_size());
		
		for(int ii = 0; ii < N; ii++){
			for(int jj = 0; jj < aa.dist().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(ii)*(aa.dist().start() + jj + 1);
			}
		}

		for(int jj = 0; jj < aa.dist().local_size(); jj++) occ[jj] = 1.0/(aa.dist().start() + jj + 1);

		auto dd = operations::calculate_density(occ, aa);
		
		for(int ii = 0; ii < M; ii++) REQUIRE(dd.linear()[ii] == Approx(0.5*ii*M*(M + 1)));
		
	}
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, M, boost::mpi3::environment::get_world_instance());

		math::array<double, 1> occ(M);
		
		for(int ii = 0; ii < N; ii++){
			for(int jj = 0; jj < aa.dist().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(ii)*(aa.dist().start() + jj + 1)*exp(complex(0.0, M_PI/65.0*ii));
			}
		}

		for(int jj = 0; jj < aa.dist().local_size(); jj++) occ[jj] = 1.0/(aa.dist().start() + jj + 1);

		auto dd = operations::calculate_density(occ, aa);
		
		for(int ii = 0; ii < M; ii++) REQUIRE(dd.linear()[ii] == Approx(0.5*ii*M*(M + 1)));
		
	}
	
}


#endif

#endif
