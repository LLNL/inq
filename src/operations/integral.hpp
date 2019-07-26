/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__INTEGRAL
#define OPERATIONS__INTEGRAL

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

#include <cassert>

namespace operations {

  template <class field_type>
  auto integral(const field_type & phi){
		//DATAOPERATIONS
		typename field_type::value_type inte = 0.0; 
		for(long ipoint = 0; ipoint < phi.basis().size(); ipoint++) inte += phi[ipoint];
		return inte*phi.basis().volume_element();
	}

  template <class field_type>
  auto diff(const field_type & phi1, const field_type & phi2){
		assert(phi1.basis() == phi2.basis());
		
		//DATAOPERATIONS
		typename field_type::value_type diff = 0.0; 
		for(long ipoint = 0; ipoint < phi1.basis().size(); ipoint++) diff += fabs(phi1[ipoint] - phi2[ipoint]);
		return diff*phi1.basis().volume_element();
	}
	
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("function operations::integral", "[integral]") {

	using namespace Catch::literals;

}


#endif
#endif
