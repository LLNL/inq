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

#include <inq_config.h>

namespace inq {
namespace math {

	template <class Type>
	class vector3 {

	public:

		vector3() = default;

		GPU_FUNCTION vector3(Type const & v0, Type const & v1, Type const & v2){
			vec_[0] = v0;
			vec_[1] = v1;
			vec_[2] = v2;
		}
		
		GPU_FUNCTION vector3(std::initializer_list<Type> const & list){
			vec_[0] = list.begin()[0];
			vec_[1] = list.begin()[1];
			vec_[2] = list.begin()[2];
		}
		
		GPU_FUNCTION vector3(std::array<Type, 3> const & arr){
			vec_[0] = arr[0];
			vec_[1] = arr[1];
			vec_[2] = arr[2];
		}
		
		GPU_FUNCTION auto & operator[](int ii){
			return vec_[ii];
		}

		GPU_FUNCTION auto & operator[](int ii) const {
			return vec_[ii];
		}

		//COMPARISON
		
		GPU_FUNCTION bool operator==(const vector3 & other) const {
			return vec_[0] == other.vec_[0] && vec_[1] == other.vec_[1] && vec_[2] == other.vec_[2];
		}

		GPU_FUNCTION bool operator!=(const vector3 & other) const {
			return ! (*this == other);
		}

		//ADDITION AND SUBSTRACTION
		
		GPU_FUNCTION vector3 & operator+=(const vector3 & other){
			vec_[0] += other.vec_[0];
			vec_[1] += other.vec_[1];
			vec_[2] += other.vec_[2];
			return *this;
		}

		GPU_FUNCTION vector3 operator+(const vector3 & other) const {
			vector3 result;
			result.vec_[0] = vec_[0] + other.vec_[0];
			result.vec_[1] = vec_[1] + other.vec_[1];
			result.vec_[2] = vec_[2] + other.vec_[2];
			return result;
		}

		GPU_FUNCTION vector3 & operator-=(const vector3 & other){
			vec_[0] -= other.vec_[0];
			vec_[1] -= other.vec_[1];
			vec_[2] -= other.vec_[2];
			return *this;
		}

		GPU_FUNCTION vector3 operator-(const vector3 & other) const {
			vector3 result;
			result.vec_[0] = vec_[0] - other.vec_[0];
			result.vec_[1] = vec_[1] - other.vec_[1];
			result.vec_[2] = vec_[2] - other.vec_[2];
			return result;
		}

		// MULTIPLICATION
		/*
		GPU_FUNCTION vector3 & operator*=(const Type & factor){
			vec_[0] *= factor;
			vec_[1] *= factor;
			vec_[2] *= factor;
			return *this;
		}
		*/
		
		//element-wise multiplication
		//		template <class TypeA, class TypeB>
		friend vector3 operator*(vector3 const & vv1, vector3 const & vv2){
			//		friend vector3<decltype(TypeA()*TypeB())> operator*(const vector3<TypeA> & vv1, const vector3<TypeB> & vv2){			
			return {vv1[0]*vv2[0], vv1[1]*vv2[1], vv1[2]*vv2[2]};
		}
		
		/*
		friend GPU_FUNCTION vector3 operator-(const vector3 & vv){
			return -1*vector3;
		}
		*/
		// INPUT OUTPUT
		
		friend std::ostream& operator <<(std::ostream & out, const vector3 & vv){
			out << vv.vec_[0] << '\t' << vv.vec_[1] << '\t' << vv.vec_[2];
			return out;
		}

		friend std::istream& operator >>(std::istream & in, vector3 & vv){
			in >> vv.vec_[0] >> vv.vec_[1] >> vv.vec_[2] ;
			return in;
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

	SECTION("Default initialization"){
		[[maybe_unused]] math::vector3<int> vv;
	}
	
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
	
	SECTION("Copy, assignment and comparison"){

		math::vector3<double> vv({0.1, 0.2, 0.3});
		math::vector3<double> vv2(vv);
		
		CHECK(vv2[0] == 0.1);
		CHECK(vv2[1] == 0.2);
		CHECK(vv2[2] == 0.3);

		math::vector3<double> vv3;

		vv3 = vv;

		CHECK(vv3[0] == 0.1);
		CHECK(vv3[1] == 0.2);
		CHECK(vv3[2] == 0.3);
		
		CHECK(vv == vv2);
		CHECK(not (vv != vv2));
		
		math::vector3<double> vv4({20.1, 0.32, 0.53});
		
		CHECK(vv4 != vv2);

		CHECK(vv4 == math::vector3<double>{20.1, 0.32, 0.53});
		
	}
		
	SECTION("Addition and substraction"){
		
		math::vector3<double> vv1({10.0, 5.0, -3.4});
		math::vector3<double> vv2({1.0, -7.8, 5.6});

		CHECK((vv1 + vv2)[0] == 11.0_a);
		CHECK((vv1 + vv2)[1] == -2.8_a);
		CHECK((vv1 + vv2)[2] ==  2.2_a);

		CHECK((vv1 - vv2)[0] ==  9.0_a);
		CHECK((vv1 - vv2)[1] == 12.8_a);
		CHECK((vv1 - vv2)[2] == -9.0_a);
		
		math::vector3<double> vv3 = vv1;

		vv3 += vv2;

		CHECK(vv3[0] == 11.0_a);
		CHECK(vv3[1] == -2.8_a);
		CHECK(vv3[2] ==  2.2_a);

		vv3 -= vv2;
		CHECK(vv3 == vv1);

	}
		
	SECTION("Multiplication"){
		
		math::vector3<double> vv1({10.0, 5.0, -3.4});
		math::vector3<double> vv2({12, -3, 4});

		auto vv3 = vv1*vv2;
		
		CHECK((vv1*vv2)[0] == 120.0_a);
		CHECK((vv1*vv2)[1] == -15.0_a);
		CHECK((vv1*vv2)[2] ==  -13.6_a);

		/*
		CHECK((vv1 - vv2)[0] ==  9.0_a);
		CHECK((vv1 - vv2)[1] == 12.8_a);
		CHECK((vv1 - vv2)[2] == -9.0_a);
		
		math::vector3<double> vv3 = vv1;

		vv3 += vv2;

		CHECK(vv3[0] == 11.0_a);
		CHECK(vv3[1] == -2.8_a);
		CHECK(vv3[2] ==  2.2_a);

		vv3 -= vv2;
		CHECK(vv3 == vv1);
		*/
	}
}

#endif
#endif
