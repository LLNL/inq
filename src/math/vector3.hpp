/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__VECTOR3
#define INQ__MATH__VECTOR3

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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
#include <math/complex.hpp>

#include <inq_config.h>

#include<array>

namespace inq {
namespace math {

	template <class Type>
	class vector3 {

	public:

		vector3() = default;

		explicit GPU_FUNCTION vector3(Type const & scal){
			vec_[0] = scal;
			vec_[1] = scal;
			vec_[2] = scal;
		}
		
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

		explicit GPU_FUNCTION vector3(Type const * const arr){
			vec_[0] = arr[0];
			vec_[1] = arr[1];
			vec_[2] = arr[2];
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

		auto data() {
			return vec_;
		}

		auto data() const {
			return vec_;
		}

		constexpr auto size() const {
			return 3;
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
		
		//element-wise multiplication and division
		template <class TypeB>
		friend GPU_FUNCTION vector3<decltype(Type()*TypeB())> operator*(const vector3 & vv1, const vector3<TypeB> & vv2){
			return {vv1[0]*vv2[0], vv1[1]*vv2[1], vv1[2]*vv2[2]};
		}

		template <class TypeB>
		friend GPU_FUNCTION vector3<decltype(TypeB()/Type())> operator/(const vector3<TypeB> & vv1, const vector3 & vv2){
			return {vv1[0]/vv2[0], vv1[1]/vv2[1], vv1[2]/vv2[2]};
		}

		//scalar multiplication and division
		private:
		template<class   > struct is_vector              : std::false_type{};
		template<class TT> struct is_vector<vector3<TT>> : std::true_type {};
		public:

		template<class TypeA, class=std::enable_if_t<not is_vector<TypeA>{}>>
		friend GPU_FUNCTION vector3<decltype(TypeA()*Type())> operator*(TypeA const& scalar, vector3 const& vv){
			return {scalar*vv[0], scalar*vv[1], scalar*vv[2]};
		}
		
		template<class TypeB, class=std::enable_if_t<not is_vector<TypeB>{}> >
		friend GPU_FUNCTION vector3<decltype(Type()*TypeB())> operator*(vector3 const& vv, TypeB const& scalar){
			return {vv[0]*scalar, vv[1]*scalar, vv[2]*scalar};
		}
		
		template <class TypeB>
		friend GPU_FUNCTION vector3<decltype(Type()/TypeB())> operator/(const vector3 & vv, const TypeB & scalar){
			return {vv[0]/scalar, vv[1]/scalar, vv[2]/scalar};
		}
		
		friend GPU_FUNCTION vector3 operator-(vector3 const & vv){
			return -1*vv;
		}

		GPU_FUNCTION vector3 & operator*=(Type const & factor){
			vec_[0] *= factor;
			vec_[1] *= factor;
			vec_[2] *= factor;
			return *this;
		}

		GPU_FUNCTION vector3 & operator/=(Type const & factor){
			vec_[0] /= factor;
			vec_[1] /= factor;
			vec_[2] /= factor;
			return *this;
		}

		// ELEMENTWISE OPERATIONS
		template <class Function>
		friend GPU_FUNCTION auto elementwise(Function const & func, vector3 const & vv) -> vector3<decltype(func(Type{}))> {
			return {func(vv[0]), func(vv[1]), func(vv[2])};
		}

		friend GPU_FUNCTION auto conj(vector3 const & vv){
			return vector3{conj(vv[0]), conj(vv[1]), conj(vv[2])};
		}

		friend GPU_FUNCTION auto real(vector3 const & vv){
			return vector3<double>{real(vv[0]), real(vv[1]), real(vv[2])};
		}
		
		friend GPU_FUNCTION auto imag(vector3 const & vv){
			return elementwise([] (auto x) { return imag(x);}, vv);
		}

		friend GPU_FUNCTION auto fabs(vector3 const & vv){
			return elementwise([] (auto x) { return fabs(x);}, vv);
		}
		
		// VECTORIAL PRODUCTS
		
		//internal product
		friend GPU_FUNCTION auto dot(vector3 const & vv1, vector3 const & vv2) {
			return conj(vv1[0])*vv2[0] + conj(vv1[1])*vv2[1] + conj(vv1[2])*vv2[2];			
		}

		//cross product
		friend GPU_FUNCTION auto cross(vector3 const & vv1, vector3 const & vv2) {
			return vector3(vv1[1]*vv2[2] - vv1[2]*vv2[1], vv1[2]*vv2[0] - vv1[0]*vv2[2], vv1[0]*vv2[1] - vv1[1]*vv2[0]);
		}
		
		friend GPU_FUNCTION auto norm(vector3 const & vv) {
			return norm(vv[0]) + norm(vv[1]) + norm(vv[2]);
		}

		friend GPU_FUNCTION auto length(vector3 const & vv) {
			return sqrt(norm(vv[0]) + norm(vv[1]) + norm(vv[2]));			
		}

		// INPUT OUTPUT
		
		friend std::ostream& operator <<(std::ostream & out, vector3 const & vv){
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

namespace std{
	template<class T> // workaround for types that recursively have a TDC workaround
	struct is_trivially_default_constructible<inq::math::vector3<T>> : 
		is_trivially_default_constructible<T>{};
}

#ifdef INQ_MATH_VECTOR3_UNIT_TEST
#undef INQ_MATH_VECTOR3_UNIT_TEST

#include <catch2/catch.hpp>

#include <iostream>

TEST_CASE("function math::vector3", "[math::vector3]") {

	using namespace inq;
	using namespace Catch::literals;

	SECTION("Default initialization"){
		math::vector3<int> vv; (void)vv;
	}

	SECTION("Scalar"){
		math::vector3<double> vv(-45.677);
		
		CHECK(vv[0] == -45.677);
		CHECK(vv[1] == -45.677);
		CHECK(vv[2] == -45.677);
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
		math::vector3<int> vv2({12, -3, 4});

		auto vv3 = vv1*vv2;
		
		CHECK(vv3[0] == 120.0_a);
		CHECK(vv3[1] == -15.0_a);
		CHECK(vv3[2] == -13.6_a);

		/*

			Disabled because std::complex causes problems on the GPU

			auto zvv = complex(0.0, 1.0)*vv1;

		CHECK(real(zvv[0]) ==   0.0_a);
		CHECK(imag(zvv[0]) ==  10.0_a);
		CHECK(real(zvv[1]) ==   0.0_a);
		CHECK(imag(zvv[1]) ==   5.0_a);
		CHECK(real(zvv[2]) ==   0.0_a);
		CHECK(imag(zvv[2]) ==  -3.4_a);

		auto zvv2 = zvv/complex(0.0, -1.0);
		zvv /= complex(0.0, -1.0);

		CHECK(real(zvv[0]) == -10.0_a);
		CHECK(imag(zvv[0]) ==   0.0_a);
		CHECK(real(zvv[1]) ==  -5.0_a);
		CHECK(imag(zvv[1]) ==   0.0_a);
		CHECK(real(zvv[2]) ==   3.4_a);
		CHECK(imag(zvv[2]) ==   0.0_a);
	
		CHECK(zvv == zvv2);
		*/
	}

	SECTION("Vector operations"){
		math::vector3<double> dv(3.0, -1.1, 0.1);

		CHECK( dot(dv, dv) == norm(dv));		
		/*

				Disabled because std::complex causes problems on the GPU

				This values might need to be updated because the conjugation was missing.

		math::vector3<complex> vv1({complex(0.0, 2.0), complex(0.2, -1.1), complex(0.1, 0.1)});
		math::vector3<complex> vv2({complex(-4.55, 9.0), complex(-0.535, -33.3), complex(2.35, -0.4)});

		CHECK((vv1|vv2) == conj(vv2|vv1));
		
		CHECK(real(vv1|vv2) == -54.462);
		CHECK(imag(vv1|vv2) == -14.9765);

		CHECK(norm(vv1) == Approx(real(conj(vv1)|vv1)));
		CHECK(imag(conj(vv1)|vv1) == 0.0_a);
		CHECK(norm(vv1) == Approx(4.0 + 0.04 + 1.21 + 0.01 + 0.01));
		CHECK(length(vv1) == Approx(sqrt(4.0 + 0.04 + 1.21 + 0.01 + 0.01)));

		CHECK(real(vv1)[0] ==  0.0_a);
		CHECK(imag(vv1)[0] ==  2.0_a);
		CHECK(real(vv1)[1] ==  0.2_a);
		CHECK(imag(vv1)[1] == -1.1_a);
		CHECK(real(vv1)[2] ==  0.1_a);
		CHECK(imag(vv1)[2] ==  0.1_a);
		*/
	}
	
	SECTION("Elementwise Multiplication same type"){

		math::vector3<double> vv1 = {10.0,  5.0, -3.4};
		math::vector3<double> vv2 = {12.0, -3.0,  4.0};

		CHECK( (vv1*vv2)[0] == vv1[0]*vv2[0] );
		CHECK( (2.*vv1)[0] == 2.*vv1[0] );
		CHECK( 2.*vv1 == vv1*2. );
	}
	
	SECTION("Elementwise Multiplication cross type"){
		/*	
				Disabled because std::complex causes problems on the GPU
				
				using complex = std::complex<double>;
				
				math::vector3<double> vv1 = {10.0,  5.0, -3.4};
				math::vector3<complex> vv2 = {12.0, -3.0,  4.0};
				
				CHECK( (vv1*vv2)[0] == vv1[0]*vv2[0] );
				CHECK( (complex{2.}*vv1)[0] == complex{2.}*vv1[0] );
		*/
	}

	SECTION("Old vector3<double> tests"){
	
		math::vector3<double> x1{1.0, 2.0, 3.0};
		math::vector3<double> x2{0.1, 0.2, 0.3};

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
	
		math::vector3<double> x3 = x1 + x2;

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
	
		CHECK(dot(x1, x2) == 1.4_a);

		auto crss = cross(x1, math::vector3<double>{-1.0, -0.5, 3.33});

		CHECK(crss[0] == 8.16_a);
		CHECK(crss[1] == -6.33_a);
		CHECK(crss[2] == 1.5_a);

		math::vector3<double> scal(6.66);

		CHECK(scal[0] == 6.66_a);
		CHECK(scal[1] == 6.66_a);
		CHECK(scal[2] == 6.66_a);

		double arr[] = {-45.0, 0.2277, 3.1};

		math::vector3<double> x4(arr);

		CHECK(x4[0] == -45.0_a);
		CHECK(x4[1] == 0.2277_a);
		CHECK(x4[2] == 3.1_a);

		x4 = math::vector3<double>(-3.3);
	
		CHECK(x4[0] == -3.3_a);
		CHECK(x4[1] == -3.3_a);
		CHECK(x4[2] == -3.3_a);

	}
	
}

#endif
#endif
