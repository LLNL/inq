/* -*- indent-tabs-mode: t -*- */

////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
//
// Produced at the Lawrence Livermore National Laboratory.
// Written by Xavier Andrade <xavier@llnl.gov>, Erik Draeger
// <draeger1@llnl.gov> and Francois Gygi <fgygi@ucdavis.edu>.
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//

//#include <config.h>
/*******************************************************************************
 *
 * spline.h
 *
 ******************************************************************************/
#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <cmath>

#define SPLINE_FLAT_BC 0.0       /* Flat boundary condition (y'=0) */
#define SPLINE_NATURAL_BC 1.e31  /* Natural boundary condition (Y"=0) */

//these functions are defined in spline.cpp, this class is an interface to them
void spline_create(const double *x, const double *y, int n, double yp1, double ypn, double *y2);
void splint (const double *xa, const double *ya, const double *y2a, int n, double x, double *y);
void splintd (const double *xa, const double *ya, const double *y2a, int n, double x, double *y, double *dy);

namespace math {
  
  class spline {
    
  public:

    enum class error {
											OUT_OF_RANGE
    };
    
    spline(){
    }
    
    spline(const double *x, double *y, int n, double yp1, double ypn){
      fit(x, y, n, yp1, ypn);
    }
    
    void fit(const double *x, double *y, int n, double yp1, double ypn){
      x_.resize(n);
      y_.resize(n);
      y2_.resize(n);
      
      for(int ii = 0; ii < n; ii++){
				x_[ii] = x[ii];
				y_[ii] = y[ii];
      }
      spline_create(x, y, n, yp1, ypn, &y2_[0]);
    }
    
    double value(const double & x) const {
      if(x < x_[0] || x > x_[x_.size() - 1]) throw error::OUT_OF_RANGE;
      
      double y;
      splint(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y);
      return y;
    }

    //OPTIMIZATION: the vectorial versions of the functions should not do the function call
    
    template <class array_type1, class array_type2>
    void value(const int size, const array_type1 & xx, array_type2 && yy) const {

      for(int ii = 0; ii < size; ii++) yy[ii] = value(xx[ii]);
    }    
    
    void derivative(const double & x, double & y, double & dy) const {
      if(x < x_[0] || x > x_[x_.size() - 1]) throw error::OUT_OF_RANGE;
      
      splintd(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y, &dy);
    }

    template <class array_type>
    void derivative(const int size, const array_type & xx, array_type & yy, array_type & dyy) const {

      for(int ii = 0; ii < size; ii++) derivative(xx[ii], yy[ii], dyy[ii]);
    }
    
    double derivative(const double & x) const {
      double y, dy;
      derivative(x, y, dy);
      return dy;
    }

    template <class array_type>
    void derivative(int size, const array_type & xx, array_type & dyy) const {
      
      for(int ii = 0; ii < size; ii++){
				double y;
				derivative(xx[ii], y, dyy[ii]);
      }
    }

    double cutoff_radius(double threshold) const {
      for(int ip = y_.size() - 1; ip >= 0; ip--){
				if(fabs(y_[ip]) >= threshold) return x_[ip];
      }
      return 0.0;
    }
    
  private :
    
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> y2_;
    
  };
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <cmath>

TEST_CASE("Class math::spline", "[spline]") {

  using namespace Catch::literals;

  const int nn = 20;
  const double dx = 0.5;
  std::vector<double> xx(nn);
  std::vector<double> ff(nn);

  for(int ii = 0; ii < nn; ii++){
    xx[ii] = dx*ii;
    ff[ii] = cos(xx[ii]);
  }

  math::spline spl;

  spl.fit(xx.data(), ff.data(), nn, SPLINE_FLAT_BC, SPLINE_NATURAL_BC);

  SECTION("Check single value interpolation"){

    REQUIRE(spl.value(0.0) == 1.0_a);
    REQUIRE(spl.value(0.4) == 0.9209862343_a);
    REQUIRE(spl.value(M_PI/2.0) == 0.0000212848_a);
    REQUIRE(spl.value(M_PI) == -0.9998862784_a);

  }

  SECTION("Check single value derivatives"){

    REQUIRE(spl.derivative(0.0) == Approx(0.0).margin(1e-8));
    REQUIRE(spl.derivative(0.4) == -0.3884027755_a);
    REQUIRE(spl.derivative(M_PI/2.0) == -0.9997597149_a);
    REQUIRE(spl.derivative(M_PI) == 0.0009543059_a);

  }
  
  SECTION("Check single value and derivative"){

    double f1, df1;

    spl.derivative(3.0/4.0*M_PI, f1, df1);
    
    REQUIRE(f1 == -0.7070426574_a);
    REQUIRE(df1 == -0.7077576301_a);

    spl.derivative(8.3, f1, df1);
    
    REQUIRE(f1 == -0.4304260777_a);
    REQUIRE(df1 == -0.9021271854_a);

  }

  SECTION("Out of range"){
    double f1, df1;
    
    REQUIRE_THROWS(spl.value(-666.6));
    REQUIRE_THROWS(spl.value(10.3));
    REQUIRE_THROWS(spl.derivative(-0.01));
    REQUIRE_THROWS(spl.derivative(1.0e10));
    REQUIRE_THROWS(spl.derivative(-0.01, f1, df1));
    REQUIRE_THROWS(spl.derivative(1000.0, f1, df1));

  }
  
  const int nn2 = 16;
  const double dx2 = 0.32;
  std::vector<double> xx2(nn2);
  std::vector<double> ff2(nn2);
  std::vector<double> dff2(nn2);
  
  //initialize the array of x values
  for(int ii = 0; ii < nn2; ii++){
    xx2[ii] = dx2*ii;
  }
  
  SECTION("Check multiple values interpolation"){

    spl.value(nn2, xx2, ff2);

    REQUIRE(ff2[13] == -0.5246734968_a);

    double diff = 0.0;
    for(int ii = 0; ii < nn2; ii++){
      diff += fabs(ff2[ii] - cos(ii*dx2));
    }
    diff /= nn2;

    REQUIRE(diff == 0.0000555057_a);
  }

  SECTION("Check multiple values derivative"){
    
    spl.derivative(nn2, xx2, dff2);
    
    REQUIRE(dff2[13] == 0.8517363378_a);

    double diff = 0.0;
    for(int ii = 0; ii < nn2; ii++){
      diff += fabs(dff2[ii] - (-sin(ii*dx2)));
    }
    diff /= nn2;

    REQUIRE(diff == 0.0004178684_a);
    
  }

  SECTION("Check multiple values interpolation and derivative"){
    
    spl.derivative(nn2, xx2, ff2, dff2);

    REQUIRE(ff2[13] == -0.5246734968_a);
    REQUIRE(dff2[13] == 0.8517363378_a);

    double diff = 0.0;
    double ddiff = 0.0;
    for(int ii = 0; ii < nn2; ii++){
      diff += fabs(ff2[ii] - cos(ii*dx2));
      ddiff += fabs(dff2[ii] - (-sin(ii*dx2)));
    }
    diff /= nn2;
    ddiff /= nn2;

    REQUIRE(diff == 0.0000555057_a);
    REQUIRE(ddiff == 0.0004178684_a);

  }

}
#endif

#endif
