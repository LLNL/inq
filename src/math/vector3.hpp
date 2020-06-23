/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__VECTOR3
#define INQ__MATH__VECTOR3

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

#ifdef HAVE_CONFIG_H
#include <inq_config.h>
#endif

namespace inq {
namespace math {

	template <class Type>
	class vector3 {

	public:

		GPU_FUNCTION vector3(){
		}
		
		GPU_FUNCTION vector3(Type const & v0, Type const & v1, Type const & v2){
			vec_[0] = v0;
			vec_[1] = v1;
			vec_[2] = v2;
		}
		
		template <class Array>
		GPU_FUNCTION vector3(Array const & arr){
			vec_[0] = arr[0];
			vec_[1] = arr[1];
			vec_[2] = arr[2];
		}
		
		GPU_FUNCTION vector3(std::initializer_list<Type> list){
			assert(list.size() == 3);
			vec_[0] = list.begin()[0];
			vec_[1] = list.begin()[1];
			vec_[2] = list.begin()[2];
		}
		
		template <class IntType>
		GPU_FUNCTION auto & operator[](IntType ii){
			return vec_[ii];
		}

		template <class IntType>
		GPU_FUNCTION auto & operator[](IntType ii) const {
			return vec_[ii];
		}
		
	private:

		Type vec_[3];

	};
 
}
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function math::vector3", "[math::vector3]") {

	using namespace inq;
	using namespace Catch::literals;

	SECTION("Initializer list"){
		math::vector3<int> vv({10, 20, 30});
		
		CHECK(vv[0] == 10);
		CHECK(vv[1] == 20);
		CHECK(vv[2] == 30);
	}

	SECTION("Initialization by element"){
		math::vector3<int> vv{1000, 2000, 3000};
		
		CHECK(vv[0] == 1000);
		CHECK(vv[1] == 2000);
		CHECK(vv[2] == 3000);
	}

	SECTION("Initialization from std::array"){

		std::array<int, 3> arr;

		arr[0] = 500;
		arr[1] = 600;
		arr[2] = 700;
		
		math::vector3<int> vv(arr);
		
		CHECK(vv[0] == 500);
		CHECK(vv[1] == 600);
		CHECK(vv[2] == 700);
	}
	
}

#endif
#endif
