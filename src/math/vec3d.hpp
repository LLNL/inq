#ifdef COMPILATION// -*-indent-tabs-mode:t;-*-
../../blds/gcc/scripts/inc++ -x c++ $0 -o $0x -lboost_serialization&&$0x&&rm $0x;exit
#endif

#ifndef VEC3D_H
#define VEC3D_H

/*
 Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa

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

#include <math/vector3.hpp>

#include <iostream>
#include <cmath>

namespace inq {
namespace math {

using vec3d = vector3<double>;

}
}

///////////////////////////////////////////////////////////////////

#if defined(INQ_UNIT_TEST) or (not __INCLUDE_LEVEL__)
#if (not __INCLUDE_LEVEL__)
#define CATCH_CONFIG_MAIN
#endif
#include <catch2/catch.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <math/array.hpp>

TEST_CASE("function math::vec3d", "[math::vec3d]") {

	using namespace inq;
	using namespace Catch::literals;
	
}
	
#endif	

#endif

