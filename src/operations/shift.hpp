/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SHIFT
#define OPERATIONS__SHIFT

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

#include <gpu/run.hpp>
#include <basis/field_set.hpp>
#include <cassert>

namespace operations {

	template <class array_1d, class field_set_type>
  void shift(const array_1d & factor, const field_set_type & shift, field_set_type & phi, double scale = 1.0){
    
    assert(size(factor) == phi.set_size());

		//DATAOPERATIONS LOOP + GPU::RUN 2D 
#ifdef HAVE_CUDA

		auto factorp = begin(factor);
		auto shiftp = begin(shift.matrix());
		auto phip = begin(phi.matrix());
		
		gpu::run(phi.set_size(), phi.basis().size(),
						 [=] __device__ (auto ist, auto ipoint){
							 phip[ipoint][ist] += scale*(factorp[ist]*shiftp[ipoint][ist]);
						 });

#else
    for(int ipoint = 0; ipoint < phi.basis().size(); ipoint++) {
			for(int ist = 0; ist < phi.set_size(); ist++) phi.matrix()[ipoint][ist] += scale*(factor[ist]*shift.matrix()[ipoint][ist]);
    }
#endif

  }
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::shift", "[operations::shift]") {

	using namespace Catch::literals;
	const int npoint = 185193;
	const int nvec = 7;
	
	basis::trivial bas(npoint);
	
	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec);
		basis::field_set<basis::trivial, double> bb(bas, nvec);

		math::array<double, 1> factor(nvec);
		
		for(int jj = 0; jj < nvec; jj++){
			for(int ii = 0; ii < npoint; ii++){
				aa.matrix()[ii][jj] = 1.0 + 0.765*ii*jj;
				bb.matrix()[ii][jj] = ii;
			}
			factor[jj] = 2.0*0.765*jj;
		}

		operations::shift(factor, bb, aa, -0.5);
				
		for(int ii = 0; ii < nvec; ii++){
			for(int jj = 0; jj < nvec; jj++) REQUIRE(aa.matrix()[ii][jj] == Approx(1.0));
		}
	}	
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec);
		basis::field_set<basis::trivial, complex> bb(bas, nvec);

		math::array<complex, 1> factor(nvec);
		
		for(int jj = 0; jj < nvec; jj++){
			for(int ii = 0; ii < npoint; ii++){
				aa.matrix()[ii][jj] = complex(ii, 1.0 + 0.765*ii*jj);
				bb.matrix()[ii][jj] = ii;
			}
			factor[jj] = complex(0.0, 2.0*0.765*jj);
		}

		operations::shift(factor, bb, aa, -0.5);
				
		for(int ii = 0; ii < npoint; ii++){
			for(int jj = 0; jj < nvec; jj++) REQUIRE(real(aa.matrix()[ii][jj]) == Approx(ii));
			for(int jj = 0; jj < nvec; jj++) REQUIRE(imag(aa.matrix()[ii][jj]) == Approx(1.0));
		}
	}	
	
	SECTION("mixed types"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec);
		basis::field_set<basis::trivial, complex> bb(bas, nvec);

		math::array<double, 1> factor(nvec);
		
		for(int jj = 0; jj < nvec; jj++){
			for(int ii = 0; ii < npoint; ii++){
				aa.matrix()[ii][jj] = complex(ii, 1.0 + 0.765*ii*jj);
				bb.matrix()[ii][jj] = complex(0.0, ii);
			}
			factor[jj] = 2.0*0.765*jj;
		}

		operations::shift(factor, bb, aa, -0.5);
				
		for(int ii = 0; ii < npoint; ii++){
			for(int jj = 0; jj < nvec; jj++) {
				REQUIRE(real(aa.matrix()[ii][jj]) == Approx(ii));
				REQUIRE(imag(aa.matrix()[ii][jj]) == Approx(1.0));
			}
		}
	}
	
}


#endif

#endif
