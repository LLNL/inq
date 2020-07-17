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


	template <class FunctionType, class VecType>
	auto finite_difference_gradient5p(FunctionType const & function, VecType const & point){
	VecType gradient = {0.0, 0.0, 0.0};
	double grid[4] = {2.0, 1.0, -1.0, -2.0};
	double coef[4] = {-1.0, 8.0, -8.0, 1.0};
	auto delta =  0.1;
	for(int idir = 0; idir < 3; idir++){
		for(int idp = 0; idp < 4; idp++ ){
			auto point_plus_delta = point;
			point_plus_delta[idir] += grid[idp]*delta;
			gradient[idir] += coef[idp]*function(point_plus_delta)/(12.0*delta);
		}
	}
	return gradient;
	}

	template <class FunctionType, class VecType>
	auto finite_difference_secondderivative5p(FunctionType const & function, VecType const & point){
	VecType divergence = {0.0, 0.0, 0.0};
	double grid[5] = { 2.0,  1.0,   0.0, -1.0, -2.0};
	double coef[5] = {-1.0, 16.0, -30.0, 16.0, -1.0};
	auto delta =  0.1;
	for(int idir = 0; idir < 3; idir++){
		for(int idp = 0; idp < 5; idp++ ){
			auto point_plus_delta = point;
			point_plus_delta[idir] += grid[idp]*delta;
			divergence[idir] += coef[idp]*function(point_plus_delta)/(12.0*delta*delta);
		}
	}
	return divergence;
	}

	template <class FunctionType, class VecType>
	auto finite_difference_divergence5p(FunctionType const & function, VecType const & point){
	auto divergence = 0.0;
	double grid[4] = {2.0, 1.0, -1.0, -2.0};
	double coef[4] = {-1.0, 8.0, -8.0, 1.0};
	auto delta =  0.1;
	for(int idir = 0; idir < 3; idir++){
		for(int idp = 0; idp < 4; idp++ ){
			auto point_plus_delta = point;
			point_plus_delta[idir] += grid[idp]*delta;
			auto f = function(point_plus_delta);
			divergence += coef[idp]*f[idir]/(12.0*delta);
		}
	}
	return divergence;
	}

}
}

///////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vec3d.hpp>

	auto gaussian_func(inq::math::vec3d rr){
		return pow(M_PI, -1.5)*exp(-norm(rr)); // sigma = 1/sqrt(2)
	}

	auto linear_func(inq::math::vec3d rr){
		return 4.5*(rr[0] + rr[1] + rr[2]);
	}

	auto quadratic_func(inq::math::vec3d rr){
		return -4.35*norm(rr);
	}

	auto dgaussian_func(inq::math::vec3d rr){
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = -2.0*rr[idir]*gaussian_func(rr);
		return ff;
	}

	auto dlinear_func(inq::math::vec3d rr){
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = 4.5;
		return ff;
	}

	auto dquadratic_func(inq::math::vec3d rr){
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = -8.7*rr[idir];
		return ff;
	}

TEST_CASE("utils::finite_difference", "[utils::finite_difference]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace inq::utils;
	using math::vec3d;

	inq::math::vec3d vec;
/*

	vec = {0.0, 0.0, 0.0};

	std::cout << finite_difference_gradient(gaussian_func, vec) << "\n";
	std::cout << finite_difference_gradient5p(gaussian_func, vec) << "  5 points \n";
	std::cout << "===============================";
	vec = {0.0001, -0.0, 0.003};

	std::cout << finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(gaussian_func, vec) - dgaussian_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(linear_func, vec) - dlinear_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(linear_func, vec) - dlinear_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
	std::cout << norm(finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec)) << "\n";
	CHECK(Approx(norm(finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec))).margin(1.0e-5) == 0.0);
	std::cout << "===============================";
	CHECK(norm(finite_difference_gradient(linear_func, vec)) == 60.75);
	CHECK(norm(finite_difference_gradient5p(linear_func, vec)) == 60.75);

	vec = {1.0, -1.0, 0.3};
	std::cout << finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(gaussian_func, vec) - dgaussian_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(linear_func, vec) - dlinear_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(linear_func, vec) - dlinear_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
//	CHECK(Approx(norm(finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec))).epsilon(0.001) == 0.0);
	CHECK(norm(finite_difference_gradient(linear_func, vec)) == 60.75);
	CHECK(norm(finite_difference_gradient5p(linear_func, vec)) == 60.75);


	vec = {0.0, 0.0, 3.3};
	std::cout << finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(gaussian_func, vec) - dgaussian_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(linear_func, vec) - dlinear_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(linear_func, vec) - dlinear_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
//	CHECK(Approx(norm(finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec))).epsilon(0.001) == 0.0);
	CHECK(norm(finite_difference_gradient(linear_func, vec)) == 60.75);
	CHECK(norm(finite_difference_gradient5p(linear_func, vec)) == 60.75);

	vec = {1000.0, -0.0020, 3.3};
	std::cout << finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(gaussian_func, vec) - dgaussian_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(linear_func, vec) - dlinear_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(linear_func, vec) - dlinear_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "\n";
	std::cout << finite_difference_gradient5p(quadratic_func, vec) - dquadratic_func(vec) << "  5 points \n";
//	CHECK(Approx(norm(finite_difference_gradient(gaussian_func, vec) - dgaussian_func(vec))).epsilon(0.001) == 0.0);
	CHECK(norm(finite_difference_gradient(linear_func, vec)) == 60.75);
	CHECK(norm(finite_difference_gradient5p(linear_func, vec)) == 60.75);
	
	*/

}

#endif
#endif



