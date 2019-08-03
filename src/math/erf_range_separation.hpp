/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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

#ifndef MATH__ERF_RANGE_SEPARATION
#define MATH__ERF_RANGE_SEPARATION

#include <cmath>

namespace math {
  class erf_range_separation {
  public:
    
    erf_range_separation(double sigma):
      sigma_(sigma){
    }

    auto long_range_potential(double rr) const {
      if(rr < 1e-8) return 2.0/(sqrt(2.0*M_PI)*sigma_);
      return erf(rr/(sigma_*sqrt(2.0)))/rr;
    }

    auto short_range_potential(double rr) const {
      return 1.0/rr - long_range_potential(rr);
    }

		[[gnu::pure]] auto short_range_potential_radius(double tol = 1e-15) const {

			// find the inverse of the function by Newton-Raphson
			double x = 1.0;
			double xold = 0.0;
			for(int ii = 0; ii < 100; ii++){
				double f = erf(x/(sigma_*sqrt(2.0))) + tol*x - 1.0;
				double df = sqrt(2.0/M_PI)*exp(-0.5*square(x/sigma_))/sigma_ + tol;
				xold = x;
				x = x - f/df;
				if(fabs(xold - x) < 1e-5) break;
			}

			return x;
    }
		
    auto long_range_density(double rr) const {
      const double exp_arg = -0.5*square(rr/sigma_);
      return exp(exp_arg)/cube(sigma_*sqrt(2.0*M_PI));
    }
    
    auto long_range_density_fourier(double gg) const {
      return exp(-0.5*square(gg*sigma_));
    }

    auto long_range_density_fourier_radius(double tol = 1e-15) const {
      return sqrt(-2.0*log(tol))/sigma_;
    }

  private:

    static double square(double x){
      return x*x;
    }

    static double cube(double x){
      return x*x*x;
    }
    
    double sigma_;
    
  };

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class math::erf_range_separation", "[math::erf_range_separation]") {

	using namespace Catch::literals;
  
  auto sep = math::erf_range_separation(0.2);

	SECTION("Potential"){

		REQUIRE(sep.short_range_potential(0.4) + sep.long_range_potential(0.4) == Approx(1.0/0.4));
		REQUIRE(sep.short_range_potential(3.3) + sep.long_range_potential(3.3) == Approx(1.0/3.3));

		REQUIRE(sep.long_range_potential(0.0) == Approx(3.989422804));
		REQUIRE(sep.long_range_potential(1.0e-7) == Approx(3.989422804));
		REQUIRE(sep.long_range_potential(1.0e-3) == Approx(3.9894061815));
		REQUIRE(sep.long_range_potential(1.0e-1) == Approx(3.8292492255));
		REQUIRE(sep.long_range_potential(5.0e-1) == Approx(1.9751613387));
		REQUIRE(sep.long_range_potential(1.0) == Approx(0.9999994267));
		REQUIRE(sep.long_range_potential(2.0) == Approx(0.5));
		REQUIRE(sep.long_range_potential(3.0) == Approx(0.3333333333));
		REQUIRE(sep.long_range_potential(4.0) == Approx(0.25));
		REQUIRE(sep.long_range_potential(5.0) == Approx(0.2));
		REQUIRE(sep.long_range_potential(6.0) == Approx(0.16666667));
		REQUIRE(sep.long_range_potential(7.0) == Approx(0.14285714));
		REQUIRE(sep.long_range_potential(8.0) == Approx(0.125));
		REQUIRE(sep.long_range_potential(9.0) == Approx(0.1111111111));
		REQUIRE(sep.long_range_potential(10.0) == Approx(0.1));

		REQUIRE(sep.short_range_potential(1.0e-7) == Approx(9999996.0105771963));
		REQUIRE(sep.short_range_potential(1.0e-3) == Approx(996.0105938185));
		REQUIRE(sep.short_range_potential(1.0e-1) == Approx(6.1707507745));
		REQUIRE(sep.short_range_potential(2.5e-1) == Approx(0.8451981893));
		REQUIRE(sep.short_range_potential(5.0e-1) == Approx(0.0248386613));
		REQUIRE(sep.short_range_potential(7.5e-1) == Approx(0.0002357794));
		REQUIRE(sep.short_range_potential(1.0) == Approx(0.0000005733));
		REQUIRE(sep.short_range_potential(1.25) == Approx(3.283621e-10));
		REQUIRE(sep.short_range_potential(1.5) == Approx(4.255855e-14));
		REQUIRE(fabs(sep.short_range_potential(1.75)) <= 1e-15);
		REQUIRE(fabs(sep.short_range_potential(2.0)) <= 1e-15);
		REQUIRE(fabs(sep.short_range_potential(2.5)) <= 1e-15);
		REQUIRE(fabs(sep.short_range_potential(3.0)) <= 1e-15);
		REQUIRE(fabs(sep.short_range_potential(3.5)) <= 1e-15);
		REQUIRE(fabs(sep.short_range_potential(4.0)) <= 1e-15);
	}

	SECTION("Short range potential radius"){

		auto radius = sep.short_range_potential_radius();
		REQUIRE(radius == 1.593078562_a);
		REQUIRE(sep.short_range_potential(radius) == 9.992007e-16_a);

		radius = sep.short_range_potential_radius(0.1);
		REQUIRE(radius == 4.089239e-01_a);
		REQUIRE(sep.short_range_potential(radius) == 0.1_a);
		
		radius = sep.short_range_potential_radius(1.0);
		REQUIRE(radius == 0.2366701558_a);
		REQUIRE(sep.short_range_potential(radius) == 1.0_a);
				
		radius = sep.short_range_potential_radius(3.0);
		REQUIRE(radius == 0.1505417002_a);
		REQUIRE(sep.short_range_potential(radius) == 3.0_a);
		
	}

	SECTION("Real space density"){

		// the integral must be 1.0
		double dx = 0.1;
		double integral = 0.0;
		for(double x = dx; x < 100.0; x += dx) integral += dx*4.0*M_PI*x*x*sep.long_range_density(x);

		REQUIRE(sep.long_range_density(0.0) == 7.9367044918_a);
		REQUIRE(integral == 1.0_a);

	}

	SECTION("Fourier space density"){

		REQUIRE(sep.long_range_density_fourier_radius() == 41.5564534067_a);
		REQUIRE(sep.long_range_density_fourier(sep.long_range_density_fourier_radius()) == 1e-15_a);

		// the G = 0 component must be 1.0
		REQUIRE(sep.long_range_density_fourier(0.0) == 1.0_a);
		
		// the integral must be zero component of the real space density
		double dg = 0.1;
		double integral = 0.0;
		for(double g = dg; g < sep.long_range_density_fourier_radius(); g += dg){
			integral += dg*4.0*M_PI*g*g*sep.long_range_density_fourier(g);
		}
		
		REQUIRE(integral/pow(2*M_PI, 3) == Approx(sep.long_range_density(0.0)));
	}
	
}

#endif

#endif

