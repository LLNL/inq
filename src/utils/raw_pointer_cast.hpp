/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__RAW_POINTER_CAST
#define INQ__UTILS__RAW_POINTER_CAST

/*
 Copyright (C) 2021 Xavier Andrade

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

namespace inq{

template <typename Type>
auto raw_pointer_cast(Type * ptr){
	return ptr;
}

}
#endif

#ifdef INQ_UTILS_RAW_POINTER_CAST_UNIT_TEST
#undef INQ_UTILS_RAW_POINTER_CAST_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <math/complex.hpp>

using namespace inq;

TEMPLATE_TEST_CASE("function raw_pointer_cast", "[raw_pointer_cast]", int, double, int const, double const, complex, complex const) {

	TestType a = 1;
	
	CHECK(typeid(raw_pointer_cast(&a)) == typeid(TestType *));
	CHECK(raw_pointer_cast(&a) == &a);
	
}
#endif

