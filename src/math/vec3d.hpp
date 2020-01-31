/* -*- indent-tabs-mode: t -*- */

#ifndef VEC3D_H
#define VEC3D_H

/*
 Copyright (C) 2020 Xavier Andrade

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cmath>
#include <cassert>

#include <gpu/run.hpp>

namespace math {
	class vec3d {
	public:

		GPU_FUNCTION vec3d(double xv, double yv, double zv) : x_{xv}, y_{yv}, z_{zv}{}
		
		GPU_FUNCTION explicit vec3d() : x_(0.0), y_(0.0), z_(0.0) {}
		
		GPU_FUNCTION explicit vec3d(const double & vv) : x_(vv), y_(vv), z_(vv) {}    

		GPU_FUNCTION explicit vec3d(const double* r) : x_(r[0]), y_(r[1]), z_(r[2]) {}

		GPU_FUNCTION double & operator[](int i){
			static_assert(sizeof(*this) == sizeof(double)*3, "must be compatible with double[3]");
			return reinterpret_cast<double*>(this)[i];
		}
		
		GPU_FUNCTION double const & operator[](int i) const{
			static_assert(sizeof(*this) == sizeof(double)*3, "must be compatible with double[3]");
			return reinterpret_cast<double const*>(this)[i];
		}
		
		GPU_FUNCTION bool operator==(const vec3d & aa) const {
			return x_ == aa.x_ && y_ == aa.y_ && z_ == aa.z_;
		}

		GPU_FUNCTION bool operator!=(const vec3d & aa) const {
			return x_ != aa.x_ || y_ != aa.y_ || z_ != aa.z_;
		}

		GPU_FUNCTION vec3d & operator +=(const vec3d & aa) {
			x_ += aa.x_;
			y_ += aa.y_;
			z_ += aa.z_;
			return *this;
		}

		GPU_FUNCTION vec3d & operator -=(const vec3d & aa) {
			x_ -= aa.x_;
			y_ -= aa.y_;
			z_ -= aa.z_;
			return *this;
		}

		GPU_FUNCTION vec3d & operator *=(const double & aa) {
			x_ *= aa;
			y_ *= aa;
			z_ *= aa;
			return *this;
		}

		GPU_FUNCTION vec3d & operator /=(const double & aa) {
			x_ /= aa;
			y_ /= aa;
			z_ /= aa;
			return *this;
		}

		GPU_FUNCTION friend const vec3d operator +(const vec3d & aa, const vec3d & bb) {
			return vec3d(aa) += bb;
		}

		GPU_FUNCTION friend const vec3d operator -(const vec3d & aa, const vec3d & bb) {
			return vec3d(aa) -= bb;
		}

		GPU_FUNCTION friend vec3d operator -(const vec3d & aa) {
			return vec3d( -aa.x_, -aa.y_, -aa.z_ );
		}

		GPU_FUNCTION friend vec3d operator *(const double & aa, const vec3d & bb) {
			return vec3d(bb) *= aa;
		}

		GPU_FUNCTION friend vec3d operator *(const vec3d & aa, const double & bb) {
			return vec3d(aa) *= bb;
		}

		GPU_FUNCTION friend vec3d operator /(const vec3d & aa, const double & bb) {
			return vec3d(aa) /= bb;
		}

		//internal product
		GPU_FUNCTION friend double operator |(const vec3d & aa, const vec3d & bb) {
			return aa.x_*bb.x_ + aa.y_*bb.y_ + aa.z_*bb.z_;
		}

		//cross product
		GPU_FUNCTION friend vec3d operator ^(const vec3d & aa, const vec3d & bb) {
			return vec3d(aa.y_*bb.z_ - aa.z_*bb.y_, aa.z_*bb.x_ - aa.x_*bb.z_, aa.x_*bb.y_ - aa.y_*bb.x_);
		}

		GPU_FUNCTION friend double norm(const vec3d& aa) {
			return aa|aa;
		}

		GPU_FUNCTION friend double length(const vec3d& aa) {
			return sqrt(aa|aa);
		}

		friend std::ostream& operator <<(std::ostream & outs, const vec3d & v) {
			outs << v.x_ << " " << v.y_ << " " << v.z_;
			return outs;
		}

		friend std::istream& operator >>(std::istream & ins, vec3d & v) {
			ins >> v.x_ >> v.y_ >> v.z_ ;
			return ins;
		}

	private:
		    
		double x_;
		double y_;
		double z_;

	};

}

///////////////////////////////////////////////////////////////////

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

#include <math/array.hpp>

TEST_CASE("function math::vec3d", "[math::vec3d]") {

	using namespace Catch::literals;

	using math::vec3d;

	vec3d x1{1.0, 2.0, 3.0};
	vec3d x2{0.1, 0.2, 0.3};
	
	REQUIRE(x1[0] == 1.0_a);
	REQUIRE(x1[1] == 2.0_a);
	REQUIRE(x1[2] == 3.0_a);

	REQUIRE(norm(x1) == 14.0_a);
	REQUIRE(length(x1) == 3.7416573868_a);
	
	REQUIRE((2.4*x1)[0] == 2.4_a);
	REQUIRE((2.4*x1)[1] == 4.8_a);
	REQUIRE((2.4*x1)[2] == 7.2_a);

	REQUIRE((x1/0.4166666666667)[0] == 2.4_a);
	REQUIRE((x1/0.4166666666667)[1] == 4.8_a);
	REQUIRE((x1/0.4166666666667)[2] == 7.2_a);

	REQUIRE(-x1 == -1.0*x1);
	
	vec3d x3 = x1 + x2;

	REQUIRE(x3[0] == 1.1_a);
	REQUIRE(x3[1] == 2.2_a);
	REQUIRE(x3[2] == 3.3_a);

	REQUIRE((x3/1.1)[0] == 1.0_a);
	REQUIRE((x3/1.1)[1] == 2.0_a);
	REQUIRE((x3/1.1)[2] == 3.0_a);

	x3 /= 2.2;

	REQUIRE(x3[0] == 0.5_a);
	REQUIRE(x3[1] == 1.0_a);
	REQUIRE(x3[2] == 1.5_a);
	
	x3 = x1 - x2;

	REQUIRE(x3[0] == 0.9_a);
	REQUIRE(x3[1] == 1.8_a);
	REQUIRE(x3[2] == 2.7_a);
	
	REQUIRE((x1|x2) == 1.4_a);

	auto cross = x1^vec3d{-1.0, -0.5, 3.33};

	REQUIRE(cross[0] == 8.16_a);
	REQUIRE(cross[1] == -6.33_a);
	REQUIRE(cross[2] == 1.5_a);

	vec3d scal(6.66);

	REQUIRE(scal[0] == 6.66_a);
	REQUIRE(scal[1] == 6.66_a);
	REQUIRE(scal[2] == 6.66_a);

	double arr[] = {-45.0, 0.2277, 3.1};

	vec3d x4(arr);

	REQUIRE(x4[0] == -45.0_a);
	REQUIRE(x4[1] == 0.2277_a);
	REQUIRE(x4[2] == 3.1_a);
	
}
	
#endif	

#endif

