/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SUM
#define OPERATIONS__SUM

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

#include <cstdlib>

namespace operations {

	template <class field_type>
	auto add(const field_type & t1, const field_type & t2){
		assert(t1.basis() == t2.basis());
		
		field_type tadd(t1.basis());
		
		//DATAOPERATIONS
		for(long ii = 0; ii < t1.basis().size(); ii++) tadd[ii] = t1[ii] + t2[ii];
		
		return tadd;
	}

	template <class field_type>
	auto add(const field_type & t1, const field_type & t2, const field_type & t3){
		assert(t1.basis() == t2.basis());
		assert(t1.basis() == t3.basis());
		
		field_type tadd(t1.basis());
		
		//DATAOPERATIONS
		for(long ii = 0; ii < t1.basis().size(); ii++) tadd[ii] = t1[ii] + t2[ii] + t3[ii];
		
		return tadd;
	}

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::add", "[add]") {

	using namespace Catch::literals;

}


#endif

#endif
