/* -*- indent-tabs-mode: t -*- */

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
#include <basis/field_set.hpp>
#include <cassert>
#include <multi/adaptors/blas.hpp>

namespace operations {

	template <class field_set_type>
  auto overlap(const field_set_type & phi1, const field_set_type & phi2){

		boost::multi::array<typename field_set_type::value_type, 2>  overlap_matrix({phi1.set_size(), phi1.set_size()});

		boost::multi::blas::gemm('N', 'C', phi1.basis().volume_element(), phi1, phi2, 0.0, overlap_matrix);

		return overlap_matrix;		
  }

	template <class field_set_type>
	auto overlap(const field_set_type & phi){

		boost::multi::array<typename field_set_type::value_type, 2>  overlap_matrix({phi.set_size(), phi.set_size()});

		// this should be done by herk, but it is not implemented yet.
		boost::multi::blas::herk('U', 'N', phi.basis().volume_element(), phi, 0.0, overlap_matrix);
	
		return overlap_matrix;
	}

	template <class field_set_type>
  auto overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

		boost::multi::array<typename field_set_type::value_type, 1>  overlap_vector(phi1.set_size());

		assert(size(overlap_vector) == phi1.set_size());

		//DATAOPERATIONS
		
		//OPTIMIZATION: this can be done more efficiently
    for(int ii = 0; ii < phi1.set_size(); ii++){
			typename field_set_type::value_type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().num_points(); ip++) aa += conj(phi1[ip][ii])*phi2[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
    }
		
		return overlap_vector;		
  }

	template <class field_set_type>
	auto overlap_diagonal(const field_set_type & phi){
		
		return overlap_diagonal(phi, phi);
	}
	
	template <class field_type>
	auto overlap_single(const field_type & phi1, const field_type & phi2){

		//DATAOPERATIONS
		//OPTIMIZATION: this can be done more efficiently
		typename field_type::value_type overlap = 0.0;
		for(int ip = 0; ip < phi1.basis().num_points(); ip++) overlap += conj(phi1[ip])*phi2[ip];
		return overlap*phi1.basis().volume_element();
	}
	
	template <class field_type>
	auto overlap_single(field_type & phi){
		return overlap_single(phi, phi);
	}
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("function operations::overlap", "[overlap]") {

	using namespace Catch::literals;

}


#endif
#endif
