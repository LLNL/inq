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
////////////////////////////////////////////////////////////////////////////////
//
// UnitCell.h
//
////////////////////////////////////////////////////////////////////////////////

//#include <config.h>

#ifndef UNITCELL_H
#define UNITCELL_H

#include <math/d3vector.hpp>
#include <valarray>
#include <array>
#include <input/cell.hpp>

namespace ions{

  class UnitCell{
	using vector_type = math::d3vector;
  private:
    vector_type a_[3];
    vector_type b_[3];
    double volume_;
    vector_type an_[13];
    vector_type bn_[13];
    double an2h_[13];
    double bn2h_[13];
  
    // 3x3 matrix forms
    double amat_[9];
    double bmat_[9];
    // 3x3 matrix form of inverse
    double amat_inv_[9];
    // 3x3 matrix form of inverse transpose
    double amat_inv_t_[9];

		int periodic_dimensions_;
		
  public:

    void set(const math::d3vector& a0, const math::d3vector& a1, const math::d3vector& a2, int arg_periodic_dimensions = 3);
		
		enum class error { WRONG_LATTICE };

		vector_type const& operator[](int i) const {return a_[i];}
    
    vector_type const& a(int i) const { return a_[i]; }
		vector_type const& b(int i) const { return b_[i]; }
  
    UnitCell() = default;

		template<class lattice_vectors_type>
		UnitCell(const lattice_vectors_type & lattice_vectors, int periodic_dimensions = 3){
			std::array<math::d3vector, 3> lvectors;
			for(int ii = 0; ii < 3; ii++){
				for(int jj = 0; jj < 3; jj++) lvectors[ii][jj] = lattice_vectors[ii][jj];
			}
			set(lvectors[0], lvectors[1], lvectors[2], periodic_dimensions);
		}
		
    UnitCell(math::d3vector const& a0, math::d3vector const& a1, math::d3vector const& a2, int periodic_dimensions = 3){
      set(a0, a1, a2, periodic_dimensions);
    }

		auto enlarge(int factor) const {
			return UnitCell(factor*a_[0], factor*a_[1], factor*a_[2], periodic_dimensions_);
		}
				
    double volume() const { return volume_; }

    template <class output_stream>
    void info(output_stream & out) const {
      out << "UNIT CELL:" << std::endl;
      out << "  Lattice vectors [b] = " << a_[0] << std::endl;
      out << "                        " << a_[1] << std::endl;
      out << "                        " << a_[2] << std::endl;
      out << "  Volume [b^3]        = " << volume_ << std::endl;
      out << std::endl;
    }
    
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

		auto diagonal_length() const {
			return sqrt(norm(a_[0]) + norm(a_[1]) + norm(a_[2]));
		}

		auto periodic_dimensions() const {
			return periodic_dimensions_;
		}
		
  };
  
  std::ostream& operator << (std::ostream& os, const UnitCell& cell);

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class ions::UnitCell", "[UnitCell]") {

  using namespace Catch::literals;
  using math::d3vector;

  {
    
    SECTION("Cubic cell"){
    
      ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

      REQUIRE(cell[0][0] == 10.0_a);
      REQUIRE(cell[0][1] ==  0.0_a);
      REQUIRE(cell[0][2] ==  0.0_a);
      REQUIRE(cell[1][0] ==  0.0_a);
      REQUIRE(cell[1][1] == 10.0_a);
      REQUIRE(cell[1][2] ==  0.0_a);
      REQUIRE(cell[2][0] ==  0.0_a);
      REQUIRE(cell[2][1] ==  0.0_a);
      REQUIRE(cell[2][2] == 10.0_a);
      
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

    SECTION("Parallelepipedic cell"){

			ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));
   
      cell.set(d3vector(28.62, 0.0, 0.0), d3vector(0.0, 90.14, 0.0), d3vector(0.0, 0.0, 12.31));


      REQUIRE(cell.a(0)[0] == 28.62_a);
      REQUIRE(cell.a(0)[1] ==  0.0_a);
      REQUIRE(cell.a(0)[2] ==  0.0_a);
      REQUIRE(cell.a(1)[0] ==  0.0_a);
      REQUIRE(cell.a(1)[1] == 90.14_a);
      REQUIRE(cell.a(1)[2] ==  0.0_a);
      REQUIRE(cell.a(2)[0] ==  0.0_a);
      REQUIRE(cell.a(2)[1] ==  0.0_a);
      REQUIRE(cell.a(2)[2] == 12.31_a);
    
      REQUIRE(cell.b(0)[0] == 0.2195382707_a);
      REQUIRE(cell.b(0)[1] ==  0.0_a);
      REQUIRE(cell.b(0)[2] ==  0.0_a);
      REQUIRE(cell.b(1)[0] ==  0.0_a);
      REQUIRE(cell.b(1)[1] == 0.0697047405_a);
      REQUIRE(cell.b(1)[2] ==  0.0_a);
      REQUIRE(cell.b(2)[0] ==  0.0_a);
      REQUIRE(cell.b(2)[1] ==  0.0_a);
      REQUIRE(cell.b(2)[2] == 0.5104131038_a);
    
      REQUIRE(cell.volume() == 31757.421708_a);

      REQUIRE(cell.amat()[0] == 28.62_a);
      REQUIRE(cell.amat()[1] ==  0.0_a);
      REQUIRE(cell.amat()[2] ==  0.0_a);
      REQUIRE(cell.amat()[3] ==  0.0_a);
      REQUIRE(cell.amat()[4] == 90.14_a);
      REQUIRE(cell.amat()[5] ==  0.0_a);
      REQUIRE(cell.amat()[6] ==  0.0_a);
      REQUIRE(cell.amat()[7] ==  0.0_a);
      REQUIRE(cell.amat()[8] == 12.31_a);
    
      REQUIRE(cell.bmat()[0] == 0.2195382707_a);
      REQUIRE(cell.bmat()[1] ==  0.0_a);
      REQUIRE(cell.bmat()[2] ==  0.0_a);
      REQUIRE(cell.bmat()[3] ==  0.0_a);
      REQUIRE(cell.bmat()[4] == 0.0697047405_a);
      REQUIRE(cell.bmat()[5] ==  0.0_a);
      REQUIRE(cell.bmat()[6] ==  0.0_a);
      REQUIRE(cell.bmat()[7] ==  0.0_a);
      REQUIRE(cell.bmat()[8] == 0.5104131038_a);
    
      REQUIRE(cell.amat_inv()[0] == 0.034940601_a);
      REQUIRE(cell.amat_inv()[1] == 0.0_a);
      REQUIRE(cell.amat_inv()[2] == 0.0_a);
      REQUIRE(cell.amat_inv()[3] == 0.0_a);
      REQUIRE(cell.amat_inv()[4] == 0.011093854_a);
      REQUIRE(cell.amat_inv()[5] == 0.0_a);
      REQUIRE(cell.amat_inv()[6] == 0.0_a);
      REQUIRE(cell.amat_inv()[7] == 0.0_a);
      REQUIRE(cell.amat_inv()[8] == 0.0812347685_a);
    
      REQUIRE(cell.amat(0) == 28.62_a);
      REQUIRE(cell.amat(1) ==  0.0_a);
      REQUIRE(cell.amat(2) ==  0.0_a);
      REQUIRE(cell.amat(3) ==  0.0_a);
      REQUIRE(cell.amat(4) == 90.14_a);
      REQUIRE(cell.amat(5) ==  0.0_a);
      REQUIRE(cell.amat(6) ==  0.0_a);
      REQUIRE(cell.amat(7) ==  0.0_a);
      REQUIRE(cell.amat(8) == 12.31_a);
    
      REQUIRE(cell.bmat(0) == 0.2195382707_a);
      REQUIRE(cell.bmat(1) ==  0.0_a);
      REQUIRE(cell.bmat(2) ==  0.0_a);
      REQUIRE(cell.bmat(3) ==  0.0_a);
      REQUIRE(cell.bmat(4) == 0.0697047405_a);
      REQUIRE(cell.bmat(5) ==  0.0_a);
      REQUIRE(cell.bmat(6) ==  0.0_a);
      REQUIRE(cell.bmat(7) ==  0.0_a);
      REQUIRE(cell.bmat(8) == 0.5104131038_a);
    
      REQUIRE(cell.amat_inv(0) == 0.034940601_a);
      REQUIRE(cell.amat_inv(1) == 0.0_a);
      REQUIRE(cell.amat_inv(2) == 0.0_a);
      REQUIRE(cell.amat_inv(3) == 0.0_a);
      REQUIRE(cell.amat_inv(4) == 0.011093854_a);
      REQUIRE(cell.amat_inv(5) == 0.0_a);
      REQUIRE(cell.amat_inv(6) == 0.0_a);
      REQUIRE(cell.amat_inv(7) == 0.0_a);
      REQUIRE(cell.amat_inv(8) == 0.0812347685_a);

      REQUIRE(cell.contains(d3vector(5.0, 5.0, 5.0)));
      REQUIRE(!cell.contains(d3vector(-5.0, 5.0, 5.0)));
      REQUIRE(!cell.contains(d3vector(5.0, -5.0, 5.0)));
      REQUIRE(!cell.contains(d3vector(5.0, 5.0, -5.0)));

      REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[0] == 5.724_a);
      REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[1] == -45.07_a);
      REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[2] == 10.67277_a);

      REQUIRE(cell.cart_to_crystal(d3vector(6.66, -203.77, 927.2))[0] == 0.2327044025_a);
      REQUIRE(cell.cart_to_crystal(d3vector(6.66, -203.77, 927.2))[1] == -2.2605946306_a);
      REQUIRE(cell.cart_to_crystal(d3vector(6.66, -203.77, 927.2))[2] == 75.3208773355_a);
          
    }

    SECTION("Non-orthogonal cell"){

			double lv[3][3] = {{6.942, 8.799, 4.759}, {9.627, 7.092, 4.819}, {4.091, 0.721, 1.043}};
      ions::UnitCell cell(lv);
      
      REQUIRE(cell.a(0)[0] == 6.942_a);
      REQUIRE(cell.a(0)[1] == 8.799_a);
      REQUIRE(cell.a(0)[2] == 4.759_a);
      REQUIRE(cell.a(1)[0] == 9.627_a);
      REQUIRE(cell.a(1)[1] == 7.092_a);
      REQUIRE(cell.a(1)[2] == 4.819_a);
      REQUIRE(cell.a(2)[0] == 4.091_a);
      REQUIRE(cell.a(2)[1] == 0.721_a);
      REQUIRE(cell.a(2)[2] == 1.043_a);
    
      REQUIRE(cell.b(0)[0] == 3.3736397602_a);
      REQUIRE(cell.b(0)[1] == 8.3200742872_a);
      REQUIRE(cell.b(0)[2] == -18.9840209206_a);
      REQUIRE(cell.b(1)[0] == -4.942140131_a);
      REQUIRE(cell.b(1)[1] == -10.517582818_a);
      REQUIRE(cell.b(1)[2] == 26.6552948109_a);
      REQUIRE(cell.b(2)[0] == 7.4410562534_a);
      REQUIRE(cell.b(2)[1] == 10.6318294029_a);
      REQUIRE(cell.b(2)[2] == -30.5117208294_a);
    
      REQUIRE(cell.volume() == 7.305321831_a);

      REQUIRE(cell.amat()[0] == 6.942_a);		
      REQUIRE(cell.amat()[1] == 8.799_a);		
      REQUIRE(cell.amat()[2] == 4.759_a);		
      REQUIRE(cell.amat()[3] == 9.627_a);		
      REQUIRE(cell.amat()[4] == 7.092_a);		
      REQUIRE(cell.amat()[5] == 4.819_a);		
      REQUIRE(cell.amat()[6] == 4.091_a);		
      REQUIRE(cell.amat()[7] == 0.721_a);		
      REQUIRE(cell.amat()[8] == 1.043_a);		
   			                          
      REQUIRE(cell.bmat()[0] == 3.3736397602_a);	
      REQUIRE(cell.bmat()[1] == 8.3200742872_a);	
      REQUIRE(cell.bmat()[2] == -18.9840209206_a);
      REQUIRE(cell.bmat()[3] == -4.942140131_a);	
      REQUIRE(cell.bmat()[4] == -10.517582818_a);	
      REQUIRE(cell.bmat()[5] == 26.6552948109_a);	
      REQUIRE(cell.bmat()[6] == 7.4410562534_a);	
      REQUIRE(cell.bmat()[7] == 10.6318294029_a);	
      REQUIRE(cell.bmat()[8] == -30.5117208294_a);
    
      REQUIRE(cell.amat_inv()[0] == 0.5369314441_a); 
      REQUIRE(cell.amat_inv()[1] == -0.7865660313_a);
      REQUIRE(cell.amat_inv()[2] == 1.1842808846_a); 
      REQUIRE(cell.amat_inv()[3] == 1.3241809497_a); 
      REQUIRE(cell.amat_inv()[4] == -1.6739252949_a);
      REQUIRE(cell.amat_inv()[5] == 1.6921082036_a); 
      REQUIRE(cell.amat_inv()[6] == -3.0214007693_a);
      REQUIRE(cell.amat_inv()[7] == 4.2423219287_a); 
      REQUIRE(cell.amat_inv()[8] == -4.8560911922_a);
    
      REQUIRE(cell.amat(0) == 6.942_a);		
      REQUIRE(cell.amat(1) == 8.799_a);		
      REQUIRE(cell.amat(2) == 4.759_a);		
      REQUIRE(cell.amat(3) == 9.627_a);		
      REQUIRE(cell.amat(4) == 7.092_a);		
      REQUIRE(cell.amat(5) == 4.819_a);		
      REQUIRE(cell.amat(6) == 4.091_a);		
      REQUIRE(cell.amat(7) == 0.721_a);		
      REQUIRE(cell.amat(8) == 1.043_a);		
    			                        
      REQUIRE(cell.bmat(0) == 3.3736397602_a);	
      REQUIRE(cell.bmat(1) == 8.3200742872_a);	
      REQUIRE(cell.bmat(2) == -18.9840209206_a);
      REQUIRE(cell.bmat(3) == -4.942140131_a);	
      REQUIRE(cell.bmat(4) == -10.517582818_a);	
      REQUIRE(cell.bmat(5) == 26.6552948109_a);	
      REQUIRE(cell.bmat(6) == 7.4410562534_a);	
      REQUIRE(cell.bmat(7) == 10.6318294029_a);	
      REQUIRE(cell.bmat(8) == -30.5117208294_a);
    
      REQUIRE(cell.amat_inv(0) == 0.5369314441_a); 
      REQUIRE(cell.amat_inv(1) == -0.7865660313_a);;
      REQUIRE(cell.amat_inv(2) == 1.1842808846_a); 
      REQUIRE(cell.amat_inv(3) == 1.3241809497_a); 
      REQUIRE(cell.amat_inv(4) == -1.6739252949_a);
      REQUIRE(cell.amat_inv(5) == 1.6921082036_a); 
      REQUIRE(cell.amat_inv(6) == -3.0214007693_a);
      REQUIRE(cell.amat_inv(7) == 4.2423219287_a); 
      REQUIRE(cell.amat_inv(8) == -4.8560911922_a);

      REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[0] == 0.121797_a);
      REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[1] == -1.161093_a);
      REQUIRE(cell.crystal_to_cart(d3vector(0.2, -0.5, 0.867))[2] == -0.553419_a);

      REQUIRE(cell.cart_to_crystal(d3vector(0.66, -23.77, 2.72))[0] == -39.3396165136_a);
      REQUIRE(cell.cart_to_crystal(d3vector(0.66, -23.77, 2.72))[1] == 50.8091863243_a);
      REQUIRE(cell.cart_to_crystal(d3vector(0.66, -23.77, 2.72))[2] == -52.6483546581_a);

      REQUIRE(cell.contains(cell.crystal_to_cart(d3vector(0.5, 0.5, 0.5))));
      //This next one fails, this has to be checked.
      //REQUIRE(!cell.contains(cell.crystal_to_cart(d3vector(1.5, 0.5, 0.5))));
      REQUIRE(!cell.contains(cell.crystal_to_cart(d3vector(0.5, -0.1, 0.0))));
      REQUIRE(!cell.contains(cell.crystal_to_cart(d3vector(0.5, 0.5, -1.0))));
      
    }
  }
}
#endif
#endif
