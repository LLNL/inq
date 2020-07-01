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

	using math::vec3d;

	vec3d x1{1.0, 2.0, 3.0};
	vec3d x2{0.1, 0.2, 0.3};

	{
		std::stringstream ss;
		boost::archive::text_oarchive{ss} << x2;
		std::cout << ss.str() <<'\n';
		vec3d x3; boost::archive::text_iarchive{ss} >> x3;
		CHECK( x3 == x2 );
	}
	
	CHECK(x1[0] == 1.0_a);
	CHECK(x1[1] == 2.0_a);
	CHECK(x1[2] == 3.0_a);

	CHECK(norm(x1) == 14.0_a);
	CHECK(length(x1) == 3.7416573868_a);
	
	CHECK((2.4*x1)[0] == 2.4_a);
	CHECK((2.4*x1)[1] == 4.8_a);
	CHECK((2.4*x1)[2] == 7.2_a);

	CHECK((x1/0.4166666666667)[0] == 2.4_a);
	CHECK((x1/0.4166666666667)[1] == 4.8_a);
	CHECK((x1/0.4166666666667)[2] == 7.2_a);

	CHECK(-x1 == -1.0*x1);
	
	vec3d x3 = x1 + x2;

	CHECK(x3[0] == 1.1_a);
	CHECK(x3[1] == 2.2_a);
	CHECK(x3[2] == 3.3_a);

	CHECK((x3/1.1)[0] == 1.0_a);
	CHECK((x3/1.1)[1] == 2.0_a);
	CHECK((x3/1.1)[2] == 3.0_a);

	x3 /= 2.2;

	CHECK(x3[0] == 0.5_a);
	CHECK(x3[1] == 1.0_a);
	CHECK(x3[2] == 1.5_a);
	
	x3 = x1 - x2;

	CHECK(x3[0] == 0.9_a);
	CHECK(x3[1] == 1.8_a);
	CHECK(x3[2] == 2.7_a);
	
	CHECK((x1|x2) == 1.4_a);

	auto cross = x1^vec3d{-1.0, -0.5, 3.33};

	CHECK(cross[0] == 8.16_a);
	CHECK(cross[1] == -6.33_a);
	CHECK(cross[2] == 1.5_a);

	vec3d scal(6.66);

	CHECK(scal[0] == 6.66_a);
	CHECK(scal[1] == 6.66_a);
	CHECK(scal[2] == 6.66_a);

	double arr[] = {-45.0, 0.2277, 3.1};

	vec3d x4(arr);

	CHECK(x4[0] == -45.0_a);
	CHECK(x4[1] == 0.2277_a);
	CHECK(x4[2] == 3.1_a);

	x4 = -3.3;
	
	CHECK(x4[0] == -3.3_a);
	CHECK(x4[1] == -3.3_a);
	CHECK(x4[2] == -3.3_a);
	
}
	
#endif	

#endif

