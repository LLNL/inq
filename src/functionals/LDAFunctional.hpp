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
////////////////////////////////////////////////////////////////////////////////
//
// LDAFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LDAFUNCTIONAL_HPP
#define LDAFUNCTIONAL_HPP

#include <vector>
#include <cassert>

namespace functionals {
  class LDAFunctional {

  public:
    
    template <class array_type1, class array_type2, class array_type3>
    static void xc_unpolarized(const int size, const array_type1 & nn, array_type2 & exc, array_type3 & vxc){
      // compute LDA xc energy and potential, unpolarized 
      // const double third=1.0/3.0;
      // c1 is (3.D0/(4.D0*pi))**third
      const double c1 = 0.6203504908994001;
      // alpha = (4/(9*pi))**third = 0.521061761198
      // const double alpha = 0.521061761198;
      // c2 = -(3/(4*pi)) / alpha = -0.458165293283
      // const double c2 = -0.458165293283;
      // c3 = (4/3) * c2 = -0.610887057711
      const double c3 = -0.610887057711;

      const double A  =  0.0311;
      const double B  = -0.048;
      const double b1 =  1.0529;
      const double b2 =  0.3334;
      const double G  = -0.1423;
  
      // C from the PZ paper: const double C  =  0.0020;
      // D from the PZ paper: const double D  = -0.0116;
      // C and D by matching Ec and Vc at rs=1
      const double D = G / ( 1.0 + b1 + b2 ) - B;
      const double C = -A - D - G * ( (b1/2.0 + b2) / ((1.0+b1+b2)*(1.0+b1+b2)));

      for(int ii = 0; ii < size; ii++){
	exc[ii] = 0.0;
	vxc[ii] = 0.0;

	if(nn[ii] > 0.0) {
	  double ro13 = cbrt(nn[ii]);
	  double rs = c1 / ro13;
	
	  double ex = 0.0, vx = 0.0, ec = 0.0, vc = 0.0;
	
	  // Next line : exchange part in Hartree units
	  vx = c3 / rs;
	  ex = 0.75 * vx;
	
	  // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
	  if ( rs < 1.0 ) {
	    double logrs = log(rs);
	    ec = A * logrs + B + C * rs * logrs + D * rs;
	    vc = A * logrs + ( B - A / 3.0 ) +
	      (2.0/3.0) * C * rs * logrs +
	      ( ( 2.0 * D - C ) / 3.0 ) * rs;
	  } else {
	    double sqrtrs = sqrt(rs);
	    double den = 1.0 + b1 * sqrtrs + b2 * rs;
	    ec = G / den;
	    vc = ec * ( 1.0 + (7.0/6.0) * b1 * sqrtrs +
			(4.0/3.0) * b2 * rs ) / den;
	  }
	  exc[ii] = ex + ec;
	  vxc[ii] = vx + vc;
	}

      }
    }

    template <class array_type1, class array_type2, class array_type3>
    static void xc_polarized(const int size, const array_type1 & nn, array_type2 & exc, array_type3 & vxc){
      // compute LDA polarized XC energy and potential

      // const double third=1.0/3.0;
      // c1 is (3.D0/(4.D0*pi))**third
      const double c1 = 0.6203504908994001;
      // alpha = (4/(9*pi))**third = 0.521061761198
      // const double alpha = 0.521061761198;
      // c2 = -(3/(4*pi)) / alpha = -0.458165293283
      // const double c2 = -0.458165293283;
      // c3 = (4/3) * c2 = -0.610887057711
      // const double c3 = -0.610887057711;
      // c4 = 2**third * c3 
      const double c4 = -0.769669463118;

      const double A =  0.01555;
      const double B = -0.0269;
      const double b1  =  1.3981;
      const double b2  =  0.2611;
      const double G   = -0.0843;
      // C from PZ paper: const double C   =  0.0007;
      // D from PZ paper: const double D   = -0.0048;
      // C and D by matching Ec and Vc at rs=1
      const double D = G / ( 1.0 + b1 + b2 ) - B;
      const double C = -A - D - G * ( (b1/2.0 + b2) / ((1.0+b1+b2)*(1.0+b1+b2)));

      for(int ii = 0; ii < size; ii++){
	exc[ii] = 0.0;
	vxc[ii] = 0.0;

	if ( nn[ii] > 0.0 ) {
	  double ro13 = cbrt(nn[ii]);
	  double rs = c1 / ro13;

	  double ex=0.0,vx=0.0,ec=0.0,vc=0.0;

	  // Next line : exchange part in Hartree units
	  vx = c4 / rs;
	  ex = 0.75 * vx;
 
	  // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
	  if ( rs < 1.0 ) {
	    double logrs = log(rs);
	    ec = A * logrs + B + C * rs * logrs + D * rs;
	    vc = A * logrs + ( B - A / 3.0 ) +
	      (2.0/3.0) * C * rs * logrs +
	      ( ( 2.0 * D - C ) / 3.0 ) * rs;
	  } else {
	    double sqrtrs = sqrt(rs);
	    double den = 1.0 + b1 * sqrtrs + b2 * rs;
	    ec = G / den;
	    vc = ec * ( 1.0 + (7.0/6.0) * b1 * sqrtrs +
			(4.0/3.0) * b2 * rs ) / den;
	  }
	  exc[ii] = ex + ec;
	  vxc[ii] = vx + vc;
	}
      }
    }
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class functionals::LDAFunctional", "[LDAFunctional]") {

  using namespace Catch::literals;

    std::vector<double> nn({1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0, 1e1, 1e2, 1e3});
    std::vector<double> exc(nn.size()), vxc(nn.size());

    //these values are from the same routine, to avoid changes in the
    //results. But they haven't been validated to other source yet.
    
    SECTION("spin unpolarized"){

      functionals::LDAFunctional::xc_unpolarized(nn.size(), nn, exc, vxc);
      
      REQUIRE(exc[0] == -0.0121328183_a);
      REQUIRE(vxc[0] == -0.0159054085_a);
      
      REQUIRE(exc[1] == -0.0246687773_a);
      REQUIRE(vxc[1] == -0.0322046187_a);
      
      REQUIRE(exc[2] == -0.0495735123_a);
      REQUIRE(vxc[2] == -0.0644958518_a);
      
      REQUIRE(exc[3] == -0.0988616346_a);
      REQUIRE(vxc[3] == -0.1284302277_a);
      
      REQUIRE(exc[4] == -0.1970983191_a);
      REQUIRE(vxc[4] == -0.2564000609_a);
      
      REQUIRE(exc[5] == -0.3962482024_a);
      REQUIRE(vxc[5] == -0.5175699500_a);
      
      REQUIRE(exc[6] == -0.8092221329_a);
      REQUIRE(vxc[6] == -1.0635879059_a);
      
      REQUIRE(exc[7] == -1.6819692861_a);
      REQUIRE(vxc[7] == -2.2215642308_a);
      
      REQUIRE(exc[8] == -3.5407734574_a);
      REQUIRE(vxc[8] == -4.6932262077_a);
      
      REQUIRE(exc[9] == -7.5211172181_a);
      REQUIRE(vxc[9] == -9.9930315852_a);
    }
    
    SECTION("spin polarized"){

      functionals::LDAFunctional::xc_polarized(nn.size(), nn, exc, vxc);

      REQUIRE(exc[0] == -0.0122936534_a);
      REQUIRE(vxc[0] == -0.0161617994_a);
      
      REQUIRE(exc[1] == -0.0253096199_a);
      REQUIRE(vxc[1] == -0.0332259764_a);
      
      REQUIRE(exc[2] == -0.0519716782_a);
      REQUIRE(vxc[2] == -0.0682116367_a);
      
      REQUIRE(exc[3] == -0.1068678070_a);
      REQUIRE(vxc[3] == -0.1404217239_a);
      
      REQUIRE(exc[4] == -0.2209158888_a);
      REQUIRE(vxc[4] == -0.2909428124_a);
      
      REQUIRE(exc[5] == -0.4603409307_a);
      REQUIRE(vxc[5] == -0.6080094243_a);
      
      REQUIRE(exc[6] == -0.968035_a);
      REQUIRE(vxc[6] == -1.28248_a);
      
      REQUIRE(exc[7] == -2.05265_a);
      REQUIRE(vxc[7] == -2.72561_a);
      
      REQUIRE(exc[8] == -4.37814_a);
      REQUIRE(vxc[8] == -5.82279_a);
      
      REQUIRE(exc[9] == -9.37581_a);
      REQUIRE(vxc[9] == -12.4826_a);
    }

}

#endif

#endif

// Local Variables:
// mode: c++
// End:
