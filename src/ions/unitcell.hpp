////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
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
////////////////////////////////////////////////////////////////////////////////
//
// UnitCell.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef UNITCELL_H
#define UNITCELL_H

#include <math/d3vector.hpp>
#include <valarray>

namespace ions{

  class UnitCell {
  private:

    math::d3vector a_[3];
    math::d3vector b_[3];
    double volume_;
    math::d3vector an_[13];
    math::d3vector bn_[13];
    double an2h_[13];
    double bn2h_[13];
  
    // 3x3 matrix forms
    double amat_[9];
    double bmat_[9];
    // 3x3 matrix form of inverse
    double amat_inv_[9];
    // 3x3 matrix form of inverse transpose
    double amat_inv_t_[9];
  
  public:

    const math::d3vector& a(int i) const { return a_[i]; }
    const math::d3vector& b(int i) const { return b_[i]; }
  
    UnitCell() { set(math::d3vector(0,0,0),math::d3vector(0,0,0),math::d3vector(0,0,0)); }  
    explicit UnitCell(const math::d3vector& a0, const math::d3vector& a1, const math::d3vector& a2)
    { set(a0,a1,a2); }
  
    void set(const math::d3vector& a0, const math::d3vector& a1, const math::d3vector& a2);
    double volume() const { return volume_; }
  
    const double* amat() const { return &amat_[0]; }
    const double* bmat() const { return &bmat_[0]; }
    const double* amat_inv() const { return &amat_inv_[0]; }
    double amat(int ij) const { return amat_[ij]; }
    double bmat(int ij) const { return bmat_[ij]; }
    double amat_inv(int ij) const { return amat_inv_[ij]; }
  
    // 3x3 matrix vector multiply Z = X Y where X is a 3x3 matrix, Y,Z 3-vectors
    void vecmult3x3(const double* x, const double* y, double *z) const;

    // 3x3 sym matrix vector multiply Z = X Y where X is a sym 3x3 matrix,
    // Y,Z 3-vectors
    void vecsmult3x3(const double* x, const double* y, double *z) const;
  
    // 3x3 matrix matrix multiply Z = X Y where X, Y are 3x3 matrices
    void matmult3x3(const double* x, const double* y, double *z) const;
    // Z = X Y where X is a symmetric 3x3 matrix and Y a general 3x3 matrix
    // uses only the first 6 elements of array xs
    // where xs[0] = x00, xs[1] = x11, xs[2] = x22,
    // xs[3] = x10, xs[4] = x21, xs[5] = x20
    void smatmult3x3(const double* xs, const double* y, double *z) const;
    void compute_deda(const std::valarray<double>& sigma, std::valarray<double>& deda) const;
  
    void cart_to_crystal(const double* scart, double* scryst) const;
    math::d3vector cart_to_crystal(const math::d3vector& v) const;
    void crystal_to_cart(const double* scryst, double* scart) const;
    math::d3vector crystal_to_cart(const math::d3vector& v) const;
    bool in_ws(const math::d3vector& v) const;
    double min_wsdist() const;
    void fold_in_ws(math::d3vector& v) const;
    bool in_bz(const math::d3vector& k) const;
    void fold_in_bz(math::d3vector& k) const;
  
    bool encloses(const UnitCell& c) const;
    bool contains(math::d3vector v) const;
  
    void print(std::ostream& os) const;  
    void printsys(std::ostream& os, std::string setcmd) const;  
    bool operator==(const UnitCell& c) const;
    bool operator!=(const UnitCell& c) const;
  };
  
  std::ostream& operator << (std::ostream& os, const UnitCell& cell);

}

#endif


#ifdef UNIT_TEST
#include <catch.hpp>

TEST_CASE("Class ions::UnitCell", "[UnitCell]") {

  using namespace Catch::literals;
  using math::d3vector;
  
  SECTION("Cubic cell"){
    
    ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

    REQUIRE(cell.a(0)[0] == 10.0_a);
    REQUIRE(cell.a(0)[1] ==  0.0_a);
    REQUIRE(cell.a(0)[2] ==  0.0_a);
    REQUIRE(cell.a(1)[0] ==  0.0_a);
    REQUIRE(cell.a(1)[1] == 10.0_a);
    REQUIRE(cell.a(1)[2] ==  0.0_a);
    REQUIRE(cell.a(2)[0] ==  0.0_a);
    REQUIRE(cell.a(2)[1] ==  0.0_a);
    REQUIRE(cell.a(2)[2] == 10.0_a);
    
    REQUIRE(cell.b(0)[0] == 0.6283185307_a);
    REQUIRE(cell.b(0)[1] ==  0.0_a);
    REQUIRE(cell.b(0)[2] ==  0.0_a);
    REQUIRE(cell.b(1)[0] ==  0.0_a);
    REQUIRE(cell.b(1)[1] == 0.6283185307_a);
    REQUIRE(cell.b(1)[2] ==  0.0_a);
    REQUIRE(cell.b(2)[0] ==  0.0_a);
    REQUIRE(cell.b(2)[1] ==  0.0_a);
    REQUIRE(cell.b(2)[2] == 0.6283185307_a);
    
    REQUIRE(cell.volume() == 1000.0_a);

    REQUIRE(cell.amat()[0] == 10.0_a);
    REQUIRE(cell.amat()[1] ==  0.0_a);
    REQUIRE(cell.amat()[2] ==  0.0_a);
    REQUIRE(cell.amat()[3] ==  0.0_a);
    REQUIRE(cell.amat()[4] == 10.0_a);
    REQUIRE(cell.amat()[5] ==  0.0_a);
    REQUIRE(cell.amat()[6] ==  0.0_a);
    REQUIRE(cell.amat()[7] ==  0.0_a);
    REQUIRE(cell.amat()[8] == 10.0_a);
    
    REQUIRE(cell.bmat()[0] == 0.6283185307_a);
    REQUIRE(cell.bmat()[1] ==  0.0_a);
    REQUIRE(cell.bmat()[2] ==  0.0_a);
    REQUIRE(cell.bmat()[3] ==  0.0_a);
    REQUIRE(cell.bmat()[4] == 0.6283185307_a);
    REQUIRE(cell.bmat()[5] ==  0.0_a);
    REQUIRE(cell.bmat()[6] ==  0.0_a);
    REQUIRE(cell.bmat()[7] ==  0.0_a);
    REQUIRE(cell.bmat()[8] == 0.6283185307_a);
    
    REQUIRE(cell.amat_inv()[0] == 0.1_a);
    REQUIRE(cell.amat_inv()[1] == 0.0_a);
    REQUIRE(cell.amat_inv()[2] == 0.0_a);
    REQUIRE(cell.amat_inv()[3] == 0.0_a);
    REQUIRE(cell.amat_inv()[4] == 0.1_a);
    REQUIRE(cell.amat_inv()[5] == 0.0_a);
    REQUIRE(cell.amat_inv()[6] == 0.0_a);
    REQUIRE(cell.amat_inv()[7] == 0.0_a);
    REQUIRE(cell.amat_inv()[8] == 0.1_a);
    
    REQUIRE(cell.amat(0) == 10.0_a);
    REQUIRE(cell.amat(1) ==  0.0_a);
    REQUIRE(cell.amat(2) ==  0.0_a);
    REQUIRE(cell.amat(3) ==  0.0_a);
    REQUIRE(cell.amat(4) == 10.0_a);
    REQUIRE(cell.amat(5) ==  0.0_a);
    REQUIRE(cell.amat(6) ==  0.0_a);
    REQUIRE(cell.amat(7) ==  0.0_a);
    REQUIRE(cell.amat(8) == 10.0_a);
    
    REQUIRE(cell.bmat(0) == 0.6283185307_a);
    REQUIRE(cell.bmat(1) ==  0.0_a);
    REQUIRE(cell.bmat(2) ==  0.0_a);
    REQUIRE(cell.bmat(3) ==  0.0_a);
    REQUIRE(cell.bmat(4) == 0.6283185307_a);
    REQUIRE(cell.bmat(5) ==  0.0_a);
    REQUIRE(cell.bmat(6) ==  0.0_a);
    REQUIRE(cell.bmat(7) ==  0.0_a);
    REQUIRE(cell.bmat(8) == 0.6283185307_a);
    
    REQUIRE(cell.amat_inv(0) == 0.1_a);
    REQUIRE(cell.amat_inv(1) == 0.0_a);
    REQUIRE(cell.amat_inv(2) == 0.0_a);
    REQUIRE(cell.amat_inv(3) == 0.0_a);
    REQUIRE(cell.amat_inv(4) == 0.1_a);
    REQUIRE(cell.amat_inv(5) == 0.0_a);
    REQUIRE(cell.amat_inv(6) == 0.0_a);
    REQUIRE(cell.amat_inv(7) == 0.0_a);
    REQUIRE(cell.amat_inv(8) == 0.1_a);

    REQUIRE(cell.contains(d3vector(5.0, 5.0, 5.0)));
    REQUIRE(!cell.contains(d3vector(-5.0, 5.0, 5.0)));
    REQUIRE(!cell.contains(d3vector(5.0, -5.0, 5.0)));
    REQUIRE(!cell.contains(d3vector(5.0, 5.0, -5.0)));

    REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[0] == 2.0_a);
    REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[1] == -5.0_a);
    REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[2] == 8.67_a);

    REQUIRE(cell.cart_to_crystal(d3vector(6.66, -3.77, 27.2))[0] == 0.666_a);
    REQUIRE(cell.cart_to_crystal(d3vector(6.66, -3.77, 27.2))[1] == -0.377_a);
    REQUIRE(cell.cart_to_crystal(d3vector(6.66, -3.77, 27.2))[2] == 2.72_a);
    
  }
  
}
#endif

// Local Variables:
// mode: c++
// End:
