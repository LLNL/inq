/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__OVERLAP
#define OPERATIONS__OVERLAP

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

#include <multi/array.hpp>
#include <states/coefficients.hpp>
#include <cassert>

namespace operations {
	
  auto overlap(const states::ks_states st, const basis::plane_wave & basis,
							 const states::coefficients & phi1, const states::coefficients & phi2) {

		boost::multi::array<states::ks_states::coeff_type, 2>  overlap_matrix({st.num_states(), st.num_states()});


		//OPTIMIZATION: this is a slow placeholder for a gemm call
    for(int ii = 0; ii < st.num_states(); ii++){
      for(int jj = 0; jj < st.num_states(); jj++){

				states::ks_states::coeff_type aa = 0.0;
				for(int kk = 0; kk < basis.num_points(); kk++) aa += std::conj(phi1.linear[ii][kk])*phi2.linear[jj][kk];
				overlap_matrix[ii][jj] = aa*basis.volume_element();

      }
    }

		return overlap_matrix;		
  }

	auto overlap(const states::ks_states st, const basis::plane_wave & basis, const states::coefficients & phi){

		//OPTIMIZATION: this can be done with syrk/herk
		return overlap(st, basis, phi, phi);
	}

	
  auto overlap_diagonal(const states::ks_states st, const basis::plane_wave & basis,
												const states::coefficients & phi1, const states::coefficients & phi2) {

		boost::multi::array<states::ks_states::coeff_type, 1>  overlap_vector(st.num_states());

		assert(size(overlap_vector) == st.num_states());
		
		//OPTIMIZATION: this can be done more efficiently
    for(int ii = 0; ii < st.num_states(); ii++){
			states::ks_states::coeff_type aa = 0.0;
			for(int kk = 0; kk < basis.num_points(); kk++) aa += std::conj(phi1.linear[ii][kk])*phi2.linear[ii][kk];
			overlap_vector[ii] = aa*basis.volume_element();
    }
		
		return overlap_vector;		
  }


	auto overlap_diagonal(const states::ks_states st, const basis::plane_wave & basis, const states::coefficients & phi){

		//OPTIMIZATION: this can be done with syrk/herk
		return overlap_diagonal(st, basis, phi, phi);
	}
	
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/plane_wave.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("function operations::overlap", "[overlap]") {

	using namespace Catch::literals;

}


#endif
#endif
