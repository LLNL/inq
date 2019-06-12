#ifndef MATH_SPHERICAL_HARMONIC
#define MATH_SPHERICAL_HARMONIC

////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov), Xavier Andrade
// (xavier@llnl.gov), and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//

namespace math {

  double spherical_harmonic(int l, int m, double gx, double gy, double gz) {
    const double pi = M_PI;
    const double fpi = 4.0 * pi;
    
    double f;
    double gnorm = sqrt(gx*gx+gy*gy+gz*gz);
    double gi = 1./gnorm;
    if (gnorm == 0.0) gi = 0.0;
    double gxi = gx*gi;
    double gyi = gy*gi;
    double gzi = gz*gi;
    double gxi2 = gxi*gxi;
    double gyi2 = gyi*gyi;
    double gzi2 = gzi*gzi;

    f = 0;
    
    // compute real spherical harmonics 
    if (l == 0) {
      const double s14pi = sqrt(1.0/fpi);
      f = s14pi;
    }
    else if (l == 1) {
      const double s34pi = sqrt(3.0/fpi);  
      if (m == 0) f = s34pi*gzi;
      else if (m == 1) f = s34pi*gxi;
      else if (m == -1) f = s34pi*gyi;
    }
    else if (l == 2) {
      const double s54pi = sqrt(5.0/fpi);
      const double s3 = sqrt(3.0);
      if (m == 0) f = 0.5*s54pi*(3.*gzi2-1.);
      else if (m == 1) f = s3*s54pi*gxi*gzi;
      else if (m == -1) f = s3*s54pi*gyi*gzi;
      else if (m == 2) f = 0.5*s3*s54pi*(gxi2-gyi2);
      else if (m == -2) f = s3*s54pi*gxi*gyi;
    }
    else if (l == 3) {
      const double s74pi = sqrt(7.0/fpi);
      const double s2132pi = sqrt(21.0/(32.*pi));
      const double s3532pi = sqrt(35.0/(32.*pi));
      const double s1054pi = sqrt(105.0/fpi);
      if (m == 0) f = 0.5*s74pi*gzi*(5.*gzi2 - 3.);
      else if (m == 1) f = s2132pi*gxi*(5.*gzi2-1.);
      else if (m == -1) f = s2132pi*gyi*(5.*gzi2-1.);
      else if (m == 2) f = 0.5*s1054pi*gzi*(gxi2-gyi2);
      else if (m == -2) f = s1054pi*gxi*gyi*gzi;
      else if (m == 3) f = s3532pi*gxi*(gxi2-3.*gyi2);
      else if (m == -3) f = s3532pi*gyi*(3.*gxi2-gyi2);
    }
    else if (l == 4) {
      const double s14pi = sqrt(1.0/fpi);
      const double s52pi = sqrt(10.0/fpi);
      const double s54pi = sqrt(5.0/fpi);
      const double s351pi = sqrt(35.0/pi);
      const double s3532pi = sqrt(35.0/(32.*pi));
      if (m == 0) f = 0.375*s14pi*(35.*gzi2*gzi2 - 30.*gzi2 + 3.);
      else if (m == 1) f = 0.75*s52pi*gxi*gzi*(7.*gzi2 - 3.);
      else if (m == -1) f = 0.75*s52pi*gyi*gzi*(7.*gzi2 - 3.);
      else if (m == 2) f = 0.1875*s54pi*(gxi2-gyi2)*(7.*gzi2 - 1.);
      else if (m == -2) f = 0.375*s54pi*gxi*gyi*(7.*gzi2 - 1.);
      else if (m == 3) f = 0.1875*s3532pi*gxi*gzi*(gxi2 - 3.*gyi2);
      else if (m == -3) f = 0.1875*s3532pi*gyi*gzi*(3.*gxi2 - gyi2);
      else if (m == 4) f = 0.1875*s351pi*gxi2*(gxi2-3.*gyi2) - gyi2*(3.*gxi2-gyi2);
      else if (m == -4) f = 0.75*s351pi*gxi*gyi*(gxi2 - gyi2);
    }
    else if (l == 5) {
      double gxi3 = gxi2*gxi;
      double gyi3 = gyi2*gyi;
      double gzi3 = gzi2*gzi;
      const double s11pi = sqrt(11.0/pi);
      const double s77pi = sqrt(77.0/pi);
      const double s1654pi = sqrt(165.0/fpi);
      const double s3852pi = sqrt(385.0/(2.*pi));
      const double s3854pi = sqrt(385.0/fpi);
      const double s11554pi = sqrt(1155.0/fpi);
      if (m == 0) f = 0.0625*s11pi*(63.*gzi3*gzi2-70.*gzi3+15.*gzi);
      else if (m == 1) f = -0.125*s1654pi*gxi*(21.*gzi2*gzi2-14.*gzi2+1.);
      else if (m == -1) f = -0.125*s1654pi*gyi*(21.*gzi2*gzi2-14.*gzi2+1.);
      else if (m == 2) f = 0.125*s11554pi*(gxi2-gyi2)*(3.*gzi3-gzi);
      else if (m == -2) f = 0.5*s11554pi*gxi*gyi*(3.*gzi3-gzi);
      else if (m == 3) f = -0.0625*s3852pi*(gxi3-3.*gxi*gyi2)*(9.*gzi2-1.);
      else if (m == -3) f = -0.0625*s3852pi*(3.*gxi2*gyi-gyi3)*(9.*gzi2-1.);
      else if (m == 4) f = 0.375*s3854pi*gzi*(gxi2*gxi2-6.*gxi2*gyi2+gyi2*gyi2);
      else if (m == -4) f = 1.5*s3854pi*gzi*(gxi3*gyi-gxi*gyi3);
      else if (m == 5) f = -0.1875*s77pi*(gxi3*gxi2-10.*gxi3*gyi2+gxi*gyi2*gyi2);
      else if (m == -5) f = -0.1875*s77pi*(5.*gxi2*gxi2*gyi-10.*gxi2*gyi3+gyi3*gyi2);
    }
    else if (l == 6) {
      double gxi3 = gxi2*gxi;
      double gyi3 = gyi2*gyi;
      double gzi3 = gzi2*gzi;
      const double s13pi = sqrt(13./pi);
      const double s2734pi = sqrt(273./fpi);
      const double s10012pi = sqrt(1001.*0.5/pi);
      const double s13652pi = sqrt(1365.*0.5/pi);
      const double s30032pi = sqrt(3003.*0.5/pi);
      const double s914pi = sqrt(91./fpi);
      if (m == 0) f = 0.03125*s13pi*(231.*gzi3*gzi3-315.*gzi2*gzi2+105.*gzi2-5.);
      else if (m == 1) f = -0.125*s2734pi*gxi*(33.*gzi3*gzi2-30.*gzi3+5.*gzi);
      else if (m == -1) f = -0.125*s2734pi*gyi*(33.*gzi3*gzi2-30.*gzi3+5.*gzi);
      else if (m == 2) f = 0.015625*s13652pi*(gxi2-gyi2)*(33.*gzi2*gzi2-18.*gzi2+1.);
      else if (m == -2) f = 0.0625*s13652pi*gxi*gyi*(33.*gzi2*gzi2-18.*gzi2+1.);
      else if (m == 3) f = -0.0625*s13652pi*(gxi3-3.*gxi*gyi2)*(11.*gzi3-3.*gzi);
      else if (m == -3) f = -0.0625*s13652pi*(3.*gxi2*gyi-gyi3)*(11.*gzi3-3.*gzi);
      else if (m == 4) f = 0.1875*s914pi*(gxi2*gxi2-6.*gxi2*gyi2+gyi2*gyi2)*(11.*gzi2-1.);
      else if (m == -4) f = 0.375*s914pi*(gxi3*gyi-gxi*gyi3)*(11.*gzi2-1.);
      else if (m == 5) f = -0.1875*s10012pi*(gxi3*gxi2-10.*gxi3*gyi2+5.*gxi*gyi2*gyi2)*gzi;
      else if (m == -5) f = -0.1875*s10012pi*(5.*gxi2*gxi2*gyi-10.*gxi2*gyi3+gyi3*gyi2)*gzi;
      else if (m == 6) f = 0.03125*s30032pi*(gxi3*gxi3-15.*gxi2*gxi2*gyi2+15.*gxi2*gyi2*gyi2-gyi3*gyi3);
      else if (m == -6) f = 0.03125*s30032pi*(6.*gxi3*gxi2*gyi-20.*gxi3*gyi3+6.*gxi*gyi3*gyi2);
    } 
    return f;
  }

  template <class array_type>
  double spherical_harmonic(int l, int m, const array_type & position){
    return spherical_harmonic(l, m, position[0], position[1], position[2]);
  }


}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <cmath>

TEST_CASE("Function math::spherical_harmonic", "[spherical_harmonic]") {

  using namespace Catch::literals;

  SECTION("l=0"){
    REQUIRE(math::spherical_harmonic(0, 0, 0.0, 0.0, 0.0) == 0.2820947918_a);

    REQUIRE(math::spherical_harmonic(0, 0, 2.0, 0.0, 0.0) == 0.2820947918_a);
    REQUIRE(math::spherical_harmonic(0, 0, 0.0, 2.0, 0.0) == 0.2820947918_a);
    REQUIRE(math::spherical_harmonic(0, 0, 0.0, 0.0, 2.0) == 0.2820947918_a);
  }

  SECTION("l=1"){
    REQUIRE(math::spherical_harmonic(1, -1, 0.0, 0.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1, -1, 2.0, 0.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1, -1, 0.0, 2.0, 0.0) == 0.4886025119_a);
    REQUIRE(math::spherical_harmonic(1, -1, 0.0, 0.0, 2.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1, -1, 0.5, 1.5, 2.5) == 0.2477666951_a);
    
    REQUIRE(math::spherical_harmonic(1,  0, 0.0, 0.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1,  0, 2.0, 0.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1,  0, 0.0, 2.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1,  0, 0.0, 0.0, 2.0) == 0.4886025119_a);
    REQUIRE(math::spherical_harmonic(1,  0, 0.5, 1.5, 2.5) == 0.4129444918_a);
    
    REQUIRE(math::spherical_harmonic(1,  1, 0.0, 0.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1,  1, 2.0, 0.0, 0.0) == 0.4886025119_a);
    REQUIRE(math::spherical_harmonic(1,  1, 0.0, 2.0, 0.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1,  1, 0.0, 0.0, 2.0) == 0.0_a);
    REQUIRE(math::spherical_harmonic(1,  1, 0.5, 1.5, 2.5) == 0.0825888984_a);
  }

}

#endif
#endif
