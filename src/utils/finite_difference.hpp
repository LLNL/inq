/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__FINITE_DIFFERENCE
#define INQ__UTILS__FINITE_DIFFERENCE

/*
 Copyright (C) 2020 Xavier Andrade, Alexey Kartsev

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

namespace inq {
namespace utils {

	// Returns value of Binomial Coefficient C(n, k)  
	int binomialCoeff(int n, int k)
	{
	if (k == 0 || k == n)  
		return 1;
	return binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k);
	}

/*
	template <class FunctionType, class VecType, class DervType >
	auto finite_difference_gradient(FunctionType const & function, VecType const & point, DervType & nk){
	double const delta = 0.01;
	VecType gradient = {0.0, 0.0, 0.0};
	for(int idir = 0; idir < 3; idir++){
		for(int in = 0; in < nk[2]+1; ix++)
		auto point_plus_delta = point;
		point_plus_delta[idir] += isign*delta;
		gradient[idir] += isign*function(point_plus_delta)/(2.0*delta);
		}
	}
	return gradient;
	}
*/

	template <class FunctionType, class VecType>
	auto finite_difference_gradient(FunctionType const & function, VecType const & point){
	double const delta = 0.01;
	VecType gradient = {0.0, 0.0, 0.0};
	for(int idir = 0; idir < 3; idir++){
		for(auto isign : {-1, 1}){
			auto point_plus_delta = point;
			point_plus_delta[idir] += isign*delta;
			gradient[idir] += isign*function(point_plus_delta)/(2.0*delta);
		}
	}
	return gradient;
	}
/*
    auto finite_difference_gradient5(FunctionType const & function, VecType const & point){
	double const delta = 0.01;
	VecType gradient = {0.0, 0.0, 0.0};
	for(int idir = 0; idir < 3; idir++){
	    for(auto isign : {-1, 1}){
		auto point_plus_delta = point;
		point_plus_delta[idir] += isign*delta;
		    gradient[idir] += isign*function(point_plus_delta)/(12.0*delta);
	    }
	}
	return gradient;
	}
*/
}
}



///////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vec3d.hpp>

	auto gaussian_func(inq::math::vec3d rr){
		return pow(M_PI, -1.5)*exp(-norm(rr)); // sigma = 1/sqrt(2)
	}

	auto dgaussian_func(inq::math::vec3d rr){
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = -2.0*rr[idir]*gaussian_func(rr);
		return ff;
	}

TEST_CASE("utils::finite_difference", "[utils::finite_difference]") {

	using namespace inq::utils;
	using namespace Catch::literals;
	using inq::math::vec3d;

	inq::math::vec3d vec;

	vec = {1, -1, 0.3};
	std::cout << finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec) << "\n";

	vec = {0, 0, 3.3};
	std::cout << finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec) << "\n";
}

#endif
#endif



