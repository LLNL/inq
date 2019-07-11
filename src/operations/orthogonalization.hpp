/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__ORTHOGONALIZATION
#define OPERATIONS__ORTHOGONALIZATION

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

#include <config.h>

#include <states/coefficients.hpp>
#include <cstdlib>

#define zpotrf FC_FUNC(zpotrf, ZPOTRF) 
extern "C" void zpotrf(const char * uplo, const int * n, complex * a, const int * lda, int * info);

#define ztrsm FC_FUNC(ztrsm, ZTRSM) 
extern "C" void ztrsm(const char * side, const char * uplo, const char * transa, const char * diag,
											const int * m, const int * n, const complex * alpha, const complex * a, const int * lda, complex * B, const int * ldb);


namespace operations {

	template <class coefficients_set_type>
  void orthogonalization(coefficients_set_type phi){

		auto olap = overlap(phi);
		
		std::cout << olap[0][0] << '\t' << olap[0][1] << '\t' << olap[0][2] << std::endl;
		std::cout << olap[1][0] << '\t' << olap[1][1] << '\t' << olap[1][2] << std::endl;
		std::cout << olap[2][0] << '\t' << olap[2][1] << '\t' << olap[2][2] << std::endl;
		
		//DATAOPERATIONS
		int info;
		const int nst = phi.set_size();
		zpotrf("U", &nst, olap.data(), &nst, &info);

		std::cout << "INFO " << info << std::endl;

		std::cout << olap[0][0] << '\t' << olap[0][1] << '\t' << olap[0][2] << std::endl;
		std::cout << olap[1][0] << '\t' << olap[1][1] << '\t' << olap[1][2] << std::endl;
		std::cout << olap[2][0] << '\t' << olap[2][1] << '\t' << olap[2][2] << std::endl;
		
		//DATAOPERATIONS
		const int np = phi.basis().num_points();
		const complex alpha = 1.0;
		ztrsm("L", "U", "T", "N", &nst, &np, &alpha, olap.data(), &nst, phi.linear.data(), &nst);
		
  }

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

#include <operations/randomize.hpp>

TEST_CASE("function operations::orthogonalization", "[orthogonalization]") {

	using namespace Catch::literals;
  using math::d3vector;

	double ecut = 25.0;
  double ll = 6.3;

	ions::geometry geo;
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::real_space pw(cell, ecut);

	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms());
	
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 6.0);

  basis::coefficients_set<basis::real_space, complex> phi(pw, 3);

	operations::randomize(phi);

	operations::orthogonalization(phi);

	auto olap = operations::overlap(phi);

	std::cout << "------" << std::endl;
	
	std::cout << olap[0][0] << '\t' << olap[0][1] << '\t' << olap[0][2] << std::endl;
	std::cout << olap[1][0] << '\t' << olap[1][1] << '\t' << olap[1][2] << std::endl;
	std::cout << olap[2][0] << '\t' << olap[2][1] << '\t' << olap[2][2] << std::endl;
	
	/*	REQUIRE(real(olap[0][0]) == 1.0_a);
	REQUIRE(imag(olap[0][0]) == 0.0_a);
	REQUIRE(real(olap[1][0]) == 0.0_a);
	REQUIRE(imag(olap[1][0]) == 0.0_a);*/
	
}


#endif

#endif
