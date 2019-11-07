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

		const auto nst = phi.set_size();
		auto occupationsp = raw_pointer_cast(occupations.data());
		auto phip = raw_pointer_cast(phi.data());
		auto densityp = raw_pointer_cast(density.data());
		
		gpu::run(phi.basis().size(),
						 [=] __device__ (auto ipoint){
							 densityp[ipoint] = 0.0;
							 for(int ist = 0; ist < nst; ist++) densityp[ipoint] += occupationsp[ist]*norm(phip[ipoint*nst + ist]);
						 });
		
#else
		
    for(int ipoint = 0; ipoint < phi.basis().size(); ipoint++){
			density[ipoint] = 0.0;
      for(int ist = 0; ist < phi.set_size(); ist++) density[ipoint] += occupations[ist]*norm(phi[ipoint][ist]);
    }
		
#endif
		
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
		
		basis::field_set<basis::trivial, double> aa(bas, M);

		math::array<double, 1> occ(M);
		
		for(int ii = 0; ii < N; ii++){
			for(int jj = 0; jj < M; jj++){
				aa[ii][jj] = sqrt(ii)*(jj + 1);
			}
		}

		for(int jj = 0; jj < M; jj++) occ[jj] = 1.0/(jj + 1);

		auto dd = operations::calculate_density(occ, aa);
		
		for(int ii = 0; ii < M; ii++) REQUIRE(dd[ii] == Approx(0.5*ii*M*(M + 1)));
		
	}
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, M);

		math::array<double, 1> occ(M);
		
		for(int ii = 0; ii < N; ii++){
			for(int jj = 0; jj < M; jj++){
				aa[ii][jj] = sqrt(ii)*(jj + 1)*exp(complex(0.0, M_PI/65.0*ii));
			}
		}

		for(int jj = 0; jj < M; jj++) occ[jj] = 1.0/(jj + 1);

		auto dd = operations::calculate_density(occ, aa);
		
		for(int ii = 0; ii < M; ii++) REQUIRE(dd[ii] == Approx(0.5*ii*M*(M + 1)));
		
	}
	
}


#endif

#endif
