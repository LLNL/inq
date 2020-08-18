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
// root directory of this partition or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// UnitCell.h
//
////////////////////////////////////////////////////////////////////////////////

//#include <inq_config.h>

#ifndef INQ__IONS__UNITCELL
#define INQ__IONS__UNITCELL

#include <math/vector3.hpp>
#include <valarray>
#include <array>
#include <input/cell.hpp>

namespace inq {
namespace ions {

  class UnitCell{
	using vector_type = math::vec3d;
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

		enum class error { WRONG_LATTICE };

    void set(const vector_type& a0, const vector_type& a1, const vector_type& a2, int arg_periodic_dimensions = 3){
			
			periodic_dimensions_ = arg_periodic_dimensions;
		
			a_[0] = a0; a_[1] = a1, a_[2] = a2;
			amat_[0] = a0[0];
			amat_[1] = a0[1];
			amat_[2] = a0[2];
			amat_[3] = a1[0];
			amat_[4] = a1[1];
			amat_[5] = a1[2];
			amat_[6] = a2[0];
			amat_[7] = a2[1];
			amat_[8] = a2[2];
  
			// volume = det(A)
			volume_ = (a0|(a1^a2));

			if(fabs(volume_) < 1e-10) throw error::WRONG_LATTICE;
		
			if ( volume_ > 0.0 ){
				// Compute rows of A-1 (columns of A^-T)
				double fac = 1.0 / volume_;
				vector_type amt0 = fac * a1 ^ a2;
				vector_type amt1 = fac * a2 ^ a0;
				vector_type amt2 = fac * a0 ^ a1;
			
				amat_inv_[0] = amt0[0];
				amat_inv_[1] = amt1[0];
				amat_inv_[2] = amt2[0];
				amat_inv_[3] = amt0[1];
				amat_inv_[4] = amt1[1];
				amat_inv_[5] = amt2[1];
				amat_inv_[6] = amt0[2];
				amat_inv_[7] = amt1[2];
				amat_inv_[8] = amt2[2];
			
				amat_inv_t_[0] = amt0[0];
				amat_inv_t_[1] = amt0[1];
				amat_inv_t_[2] = amt0[2];
				amat_inv_t_[3] = amt1[0];
				amat_inv_t_[4] = amt1[1];
				amat_inv_t_[5] = amt1[2];
				amat_inv_t_[6] = amt2[0];
				amat_inv_t_[7] = amt2[1];
				amat_inv_t_[8] = amt2[2];
			
				// B = 2 pi A^-T
				b_[0] = 2.0 * M_PI * amt0;
				b_[1] = 2.0 * M_PI * amt1;
				b_[2] = 2.0 * M_PI * amt2;
			
				bmat_[0] = b_[0][0];
				bmat_[1] = b_[0][1];
				bmat_[2] = b_[0][2];
				bmat_[3] = b_[1][0];
				bmat_[4] = b_[1][1];
				bmat_[5] = b_[1][2];
				bmat_[6] = b_[2][0];
				bmat_[7] = b_[2][1];
				bmat_[8] = b_[2][2];
			} else  {
				b_[0] = b_[1] = b_[2] = vector_type(0.0,0.0,0.0);
				amat_inv_[0] =  amat_inv_[1] =  amat_inv_[2] = 
					amat_inv_[3] =  amat_inv_[4] =  amat_inv_[5] = 
					amat_inv_[6] =  amat_inv_[7] =  amat_inv_[8] = 0.0;
				bmat_[0] =  bmat_[1] =  bmat_[2] = 
					bmat_[3] =  bmat_[4] =  bmat_[5] = 
					bmat_[6] =  bmat_[7] =  bmat_[8] = 0.0;
			}
		
			an_[0]  = a_[0];
			an_[1]  = a_[1];
			an_[2]  = a_[2];
			an_[3]  = a_[0] + a_[1];
			an_[4]  = a_[0] - a_[1];
			an_[5]  = a_[1] + a_[2];
			an_[6]  = a_[1] - a_[2];
			an_[7]  = a_[2] + a_[0];
			an_[8]  = a_[2] - a_[0];
			an_[9]  = a_[0] + a_[1] + a_[2];
			an_[10] = a_[0] - a_[1] - a_[2];
			an_[11] = a_[0] + a_[1] - a_[2];
			an_[12] = a_[0] - a_[1] + a_[2];
 
			bn_[0]  = b_[0];
			bn_[1]  = b_[1];
			bn_[2]  = b_[2];
			bn_[3]  = b_[0] + b_[1];
			bn_[4]  = b_[0] - b_[1];
			bn_[5]  = b_[1] + b_[2];
			bn_[6]  = b_[1] - b_[2];
			bn_[7]  = b_[2] + b_[0];
			bn_[8]  = b_[2] - b_[0];
			bn_[9]  = b_[0] + b_[1] + b_[2];
			bn_[10] = b_[0] - b_[1] - b_[2];
			bn_[11] = b_[0] + b_[1] - b_[2];
			bn_[12] = b_[0] - b_[1] + b_[2];
 
			for(int i = 0; i < 13; i++ ){
				an2h_[i] = 0.5 * norm(an_[i]);
				bn2h_[i] = 0.5 * norm(bn_[i]);
			}

		}
		
		vector_type const& operator[](int i) const {return a_[i];}
    
    vector_type const& a(int i) const { return a_[i]; }
		vector_type const& b(int i) const { return b_[i]; }
  
    UnitCell() = default;

		template<class lattice_vectors_type>
		UnitCell(const lattice_vectors_type & lattice_vectors, int periodic_dimensions = 3){
			std::array<math::vec3d, 3> lvectors;
			for(int ii = 0; ii < 3; ii++){
				for(int jj = 0; jj < 3; jj++) lvectors[ii][jj] = lattice_vectors[ii][jj];
			}
			set(lvectors[0], lvectors[1], lvectors[2], periodic_dimensions);
		}
		
    UnitCell(math::vec3d const& a0, math::vec3d const& a1, math::vec3d const& a2, int periodic_dimensions = 3){
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

    math::vec3d cart_to_crystal(const math::vec3d& v) const {
			vector_type vcryst;
			const double twopiinv = 0.5/M_PI;
			vcryst[0] = (b_[0]|v)*twopiinv;
			vcryst[1] = (b_[1]|v)*twopiinv;
			vcryst[2] = (b_[2]|v)*twopiinv;
			return vcryst;
		}

    void crystal_to_cart(const double* scryst, double* scart) const;

    math::vec3d crystal_to_cart(const math::vec3d& v) const {
			vec3d vcart = v[0]*a_[0] + v[1]*a_[1] + v[2]*a_[2];
			return vcart;
		}
		
    bool in_ws(const math::vec3d& v) const;
    double min_wsdist() const;
    void fold_in_ws(math::vec3d& v) const;
    bool in_bz(const math::vec3d& k) const;
    void fold_in_bz(math::vec3d& k) const;
  
    bool encloses(const UnitCell& c) const;
    bool contains(math::vec3d v) const;
  
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
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class ions::UnitCell", "[UnitCell]") {

	using namespace inq;
	using namespace Catch::literals;
  using math::vec3d;

  {
    
    SECTION("Cubic cell"){
    
      ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 10.0, 0.0), vec3d(0.0, 0.0, 10.0));

      CHECK(cell[0][0] == 10.0_a);
      CHECK(cell[0][1] ==  0.0_a);
      CHECK(cell[0][2] ==  0.0_a);
      CHECK(cell[1][0] ==  0.0_a);
      CHECK(cell[1][1] == 10.0_a);
      CHECK(cell[1][2] ==  0.0_a);
      CHECK(cell[2][0] ==  0.0_a);
      CHECK(cell[2][1] ==  0.0_a);
      CHECK(cell[2][2] == 10.0_a);
      
      CHECK(cell.a(0)[0] == 10.0_a);
      CHECK(cell.a(0)[1] ==  0.0_a);
      CHECK(cell.a(0)[2] ==  0.0_a);
      CHECK(cell.a(1)[0] ==  0.0_a);
      CHECK(cell.a(1)[1] == 10.0_a);
      CHECK(cell.a(1)[2] ==  0.0_a);
      CHECK(cell.a(2)[0] ==  0.0_a);
      CHECK(cell.a(2)[1] ==  0.0_a);
      CHECK(cell.a(2)[2] == 10.0_a);
    
      CHECK(cell.b(0)[0] == 0.6283185307_a);
      CHECK(cell.b(0)[1] ==  0.0_a);
      CHECK(cell.b(0)[2] ==  0.0_a);
      CHECK(cell.b(1)[0] ==  0.0_a);
      CHECK(cell.b(1)[1] == 0.6283185307_a);
      CHECK(cell.b(1)[2] ==  0.0_a);
      CHECK(cell.b(2)[0] ==  0.0_a);
      CHECK(cell.b(2)[1] ==  0.0_a);
      CHECK(cell.b(2)[2] == 0.6283185307_a);
    
      CHECK(cell.volume() == 1000.0_a);

      CHECK(cell.amat()[0] == 10.0_a);
      CHECK(cell.amat()[1] ==  0.0_a);
      CHECK(cell.amat()[2] ==  0.0_a);
      CHECK(cell.amat()[3] ==  0.0_a);
      CHECK(cell.amat()[4] == 10.0_a);
      CHECK(cell.amat()[5] ==  0.0_a);
      CHECK(cell.amat()[6] ==  0.0_a);
      CHECK(cell.amat()[7] ==  0.0_a);
      CHECK(cell.amat()[8] == 10.0_a);
    
      CHECK(cell.bmat()[0] == 0.6283185307_a);
      CHECK(cell.bmat()[1] ==  0.0_a);
      CHECK(cell.bmat()[2] ==  0.0_a);
      CHECK(cell.bmat()[3] ==  0.0_a);
      CHECK(cell.bmat()[4] == 0.6283185307_a);
      CHECK(cell.bmat()[5] ==  0.0_a);
      CHECK(cell.bmat()[6] ==  0.0_a);
      CHECK(cell.bmat()[7] ==  0.0_a);
      CHECK(cell.bmat()[8] == 0.6283185307_a);
    
      CHECK(cell.amat_inv()[0] == 0.1_a);
      CHECK(cell.amat_inv()[1] == 0.0_a);
      CHECK(cell.amat_inv()[2] == 0.0_a);
      CHECK(cell.amat_inv()[3] == 0.0_a);
      CHECK(cell.amat_inv()[4] == 0.1_a);
      CHECK(cell.amat_inv()[5] == 0.0_a);
      CHECK(cell.amat_inv()[6] == 0.0_a);
      CHECK(cell.amat_inv()[7] == 0.0_a);
      CHECK(cell.amat_inv()[8] == 0.1_a);
    
      CHECK(cell.amat(0) == 10.0_a);
      CHECK(cell.amat(1) ==  0.0_a);
      CHECK(cell.amat(2) ==  0.0_a);
      CHECK(cell.amat(3) ==  0.0_a);
      CHECK(cell.amat(4) == 10.0_a);
      CHECK(cell.amat(5) ==  0.0_a);
      CHECK(cell.amat(6) ==  0.0_a);
      CHECK(cell.amat(7) ==  0.0_a);
      CHECK(cell.amat(8) == 10.0_a);
    
      CHECK(cell.bmat(0) == 0.6283185307_a);
      CHECK(cell.bmat(1) ==  0.0_a);
      CHECK(cell.bmat(2) ==  0.0_a);
      CHECK(cell.bmat(3) ==  0.0_a);
      CHECK(cell.bmat(4) == 0.6283185307_a);
      CHECK(cell.bmat(5) ==  0.0_a);
      CHECK(cell.bmat(6) ==  0.0_a);
      CHECK(cell.bmat(7) ==  0.0_a);
      CHECK(cell.bmat(8) == 0.6283185307_a);
    
      CHECK(cell.amat_inv(0) == 0.1_a);
      CHECK(cell.amat_inv(1) == 0.0_a);
      CHECK(cell.amat_inv(2) == 0.0_a);
      CHECK(cell.amat_inv(3) == 0.0_a);
      CHECK(cell.amat_inv(4) == 0.1_a);
      CHECK(cell.amat_inv(5) == 0.0_a);
      CHECK(cell.amat_inv(6) == 0.0_a);
      CHECK(cell.amat_inv(7) == 0.0_a);
      CHECK(cell.amat_inv(8) == 0.1_a);

      CHECK(cell.contains(vec3d(5.0, 5.0, 5.0)));
      CHECK(!cell.contains(vec3d(-5.0, 5.0, 5.0)));
      CHECK(!cell.contains(vec3d(5.0, -5.0, 5.0)));
      CHECK(!cell.contains(vec3d(5.0, 5.0, -5.0)));

      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[0] == 2.0_a);
      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[1] == -5.0_a);
      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[2] == 8.67_a);

      CHECK(cell.cart_to_crystal(vec3d(6.66, -3.77, 27.2))[0] == 0.666_a);
      CHECK(cell.cart_to_crystal(vec3d(6.66, -3.77, 27.2))[1] == -0.377_a);
      CHECK(cell.cart_to_crystal(vec3d(6.66, -3.77, 27.2))[2] == 2.72_a);
    
    }

    SECTION("Parallelepipedic cell"){

			ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 10.0, 0.0), vec3d(0.0, 0.0, 10.0));
   
      cell.set(vec3d(28.62, 0.0, 0.0), vec3d(0.0, 90.14, 0.0), vec3d(0.0, 0.0, 12.31));


      CHECK(cell.a(0)[0] == 28.62_a);
      CHECK(cell.a(0)[1] ==  0.0_a);
      CHECK(cell.a(0)[2] ==  0.0_a);
      CHECK(cell.a(1)[0] ==  0.0_a);
      CHECK(cell.a(1)[1] == 90.14_a);
      CHECK(cell.a(1)[2] ==  0.0_a);
      CHECK(cell.a(2)[0] ==  0.0_a);
      CHECK(cell.a(2)[1] ==  0.0_a);
      CHECK(cell.a(2)[2] == 12.31_a);
    
      CHECK(cell.b(0)[0] == 0.2195382707_a);
      CHECK(cell.b(0)[1] ==  0.0_a);
      CHECK(cell.b(0)[2] ==  0.0_a);
      CHECK(cell.b(1)[0] ==  0.0_a);
      CHECK(cell.b(1)[1] == 0.0697047405_a);
      CHECK(cell.b(1)[2] ==  0.0_a);
      CHECK(cell.b(2)[0] ==  0.0_a);
      CHECK(cell.b(2)[1] ==  0.0_a);
      CHECK(cell.b(2)[2] == 0.5104131038_a);
    
      CHECK(cell.volume() == 31757.421708_a);

      CHECK(cell.amat()[0] == 28.62_a);
      CHECK(cell.amat()[1] ==  0.0_a);
      CHECK(cell.amat()[2] ==  0.0_a);
      CHECK(cell.amat()[3] ==  0.0_a);
      CHECK(cell.amat()[4] == 90.14_a);
      CHECK(cell.amat()[5] ==  0.0_a);
      CHECK(cell.amat()[6] ==  0.0_a);
      CHECK(cell.amat()[7] ==  0.0_a);
      CHECK(cell.amat()[8] == 12.31_a);
    
      CHECK(cell.bmat()[0] == 0.2195382707_a);
      CHECK(cell.bmat()[1] ==  0.0_a);
      CHECK(cell.bmat()[2] ==  0.0_a);
      CHECK(cell.bmat()[3] ==  0.0_a);
      CHECK(cell.bmat()[4] == 0.0697047405_a);
      CHECK(cell.bmat()[5] ==  0.0_a);
      CHECK(cell.bmat()[6] ==  0.0_a);
      CHECK(cell.bmat()[7] ==  0.0_a);
      CHECK(cell.bmat()[8] == 0.5104131038_a);
    
      CHECK(cell.amat_inv()[0] == 0.034940601_a);
      CHECK(cell.amat_inv()[1] == 0.0_a);
      CHECK(cell.amat_inv()[2] == 0.0_a);
      CHECK(cell.amat_inv()[3] == 0.0_a);
      CHECK(cell.amat_inv()[4] == 0.011093854_a);
      CHECK(cell.amat_inv()[5] == 0.0_a);
      CHECK(cell.amat_inv()[6] == 0.0_a);
      CHECK(cell.amat_inv()[7] == 0.0_a);
      CHECK(cell.amat_inv()[8] == 0.0812347685_a);
    
      CHECK(cell.amat(0) == 28.62_a);
      CHECK(cell.amat(1) ==  0.0_a);
      CHECK(cell.amat(2) ==  0.0_a);
      CHECK(cell.amat(3) ==  0.0_a);
      CHECK(cell.amat(4) == 90.14_a);
      CHECK(cell.amat(5) ==  0.0_a);
      CHECK(cell.amat(6) ==  0.0_a);
      CHECK(cell.amat(7) ==  0.0_a);
      CHECK(cell.amat(8) == 12.31_a);
    
      CHECK(cell.bmat(0) == 0.2195382707_a);
      CHECK(cell.bmat(1) ==  0.0_a);
      CHECK(cell.bmat(2) ==  0.0_a);
      CHECK(cell.bmat(3) ==  0.0_a);
      CHECK(cell.bmat(4) == 0.0697047405_a);
      CHECK(cell.bmat(5) ==  0.0_a);
      CHECK(cell.bmat(6) ==  0.0_a);
      CHECK(cell.bmat(7) ==  0.0_a);
      CHECK(cell.bmat(8) == 0.5104131038_a);
    
      CHECK(cell.amat_inv(0) == 0.034940601_a);
      CHECK(cell.amat_inv(1) == 0.0_a);
      CHECK(cell.amat_inv(2) == 0.0_a);
      CHECK(cell.amat_inv(3) == 0.0_a);
      CHECK(cell.amat_inv(4) == 0.011093854_a);
      CHECK(cell.amat_inv(5) == 0.0_a);
      CHECK(cell.amat_inv(6) == 0.0_a);
      CHECK(cell.amat_inv(7) == 0.0_a);
      CHECK(cell.amat_inv(8) == 0.0812347685_a);

      CHECK(cell.contains(vec3d(5.0, 5.0, 5.0)));
      CHECK(!cell.contains(vec3d(-5.0, 5.0, 5.0)));
      CHECK(!cell.contains(vec3d(5.0, -5.0, 5.0)));
      CHECK(!cell.contains(vec3d(5.0, 5.0, -5.0)));

      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[0] == 5.724_a);
      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[1] == -45.07_a);
      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[2] == 10.67277_a);

      CHECK(cell.cart_to_crystal(vec3d(6.66, -203.77, 927.2))[0] == 0.2327044025_a);
      CHECK(cell.cart_to_crystal(vec3d(6.66, -203.77, 927.2))[1] == -2.2605946306_a);
      CHECK(cell.cart_to_crystal(vec3d(6.66, -203.77, 927.2))[2] == 75.3208773355_a);
          
    }

    SECTION("Non-orthogonal cell"){

			double lv[3][3] = {{6.942, 8.799, 4.759}, {9.627, 7.092, 4.819}, {4.091, 0.721, 1.043}};
      ions::UnitCell cell(lv);
      
      CHECK(cell.a(0)[0] == 6.942_a);
      CHECK(cell.a(0)[1] == 8.799_a);
      CHECK(cell.a(0)[2] == 4.759_a);
      CHECK(cell.a(1)[0] == 9.627_a);
      CHECK(cell.a(1)[1] == 7.092_a);
      CHECK(cell.a(1)[2] == 4.819_a);
      CHECK(cell.a(2)[0] == 4.091_a);
      CHECK(cell.a(2)[1] == 0.721_a);
      CHECK(cell.a(2)[2] == 1.043_a);
    
      CHECK(cell.b(0)[0] == 3.3736397602_a);
      CHECK(cell.b(0)[1] == 8.3200742872_a);
      CHECK(cell.b(0)[2] == -18.9840209206_a);
      CHECK(cell.b(1)[0] == -4.942140131_a);
      CHECK(cell.b(1)[1] == -10.517582818_a);
      CHECK(cell.b(1)[2] == 26.6552948109_a);
      CHECK(cell.b(2)[0] == 7.4410562534_a);
      CHECK(cell.b(2)[1] == 10.6318294029_a);
      CHECK(cell.b(2)[2] == -30.5117208294_a);
    
      CHECK(cell.volume() == 7.305321831_a);

      CHECK(cell.amat()[0] == 6.942_a);		
      CHECK(cell.amat()[1] == 8.799_a);		
      CHECK(cell.amat()[2] == 4.759_a);		
      CHECK(cell.amat()[3] == 9.627_a);		
      CHECK(cell.amat()[4] == 7.092_a);		
      CHECK(cell.amat()[5] == 4.819_a);		
      CHECK(cell.amat()[6] == 4.091_a);		
      CHECK(cell.amat()[7] == 0.721_a);		
      CHECK(cell.amat()[8] == 1.043_a);		
   			                          
      CHECK(cell.bmat()[0] == 3.3736397602_a);	
      CHECK(cell.bmat()[1] == 8.3200742872_a);	
      CHECK(cell.bmat()[2] == -18.9840209206_a);
      CHECK(cell.bmat()[3] == -4.942140131_a);	
      CHECK(cell.bmat()[4] == -10.517582818_a);	
      CHECK(cell.bmat()[5] == 26.6552948109_a);	
      CHECK(cell.bmat()[6] == 7.4410562534_a);	
      CHECK(cell.bmat()[7] == 10.6318294029_a);	
      CHECK(cell.bmat()[8] == -30.5117208294_a);
    
      CHECK(cell.amat_inv()[0] == 0.5369314441_a); 
      CHECK(cell.amat_inv()[1] == -0.7865660313_a);
      CHECK(cell.amat_inv()[2] == 1.1842808846_a); 
      CHECK(cell.amat_inv()[3] == 1.3241809497_a); 
      CHECK(cell.amat_inv()[4] == -1.6739252949_a);
      CHECK(cell.amat_inv()[5] == 1.6921082036_a); 
      CHECK(cell.amat_inv()[6] == -3.0214007693_a);
      CHECK(cell.amat_inv()[7] == 4.2423219287_a); 
      CHECK(cell.amat_inv()[8] == -4.8560911922_a);
    
      CHECK(cell.amat(0) == 6.942_a);		
      CHECK(cell.amat(1) == 8.799_a);		
      CHECK(cell.amat(2) == 4.759_a);		
      CHECK(cell.amat(3) == 9.627_a);		
      CHECK(cell.amat(4) == 7.092_a);		
      CHECK(cell.amat(5) == 4.819_a);		
      CHECK(cell.amat(6) == 4.091_a);		
      CHECK(cell.amat(7) == 0.721_a);		
      CHECK(cell.amat(8) == 1.043_a);		
    			                        
      CHECK(cell.bmat(0) == 3.3736397602_a);	
      CHECK(cell.bmat(1) == 8.3200742872_a);	
      CHECK(cell.bmat(2) == -18.9840209206_a);
      CHECK(cell.bmat(3) == -4.942140131_a);	
      CHECK(cell.bmat(4) == -10.517582818_a);	
      CHECK(cell.bmat(5) == 26.6552948109_a);	
      CHECK(cell.bmat(6) == 7.4410562534_a);	
      CHECK(cell.bmat(7) == 10.6318294029_a);	
      CHECK(cell.bmat(8) == -30.5117208294_a);
    
      CHECK(cell.amat_inv(0) == 0.5369314441_a); 
      CHECK(cell.amat_inv(1) == -0.7865660313_a);;
      CHECK(cell.amat_inv(2) == 1.1842808846_a); 
      CHECK(cell.amat_inv(3) == 1.3241809497_a); 
      CHECK(cell.amat_inv(4) == -1.6739252949_a);
      CHECK(cell.amat_inv(5) == 1.6921082036_a); 
      CHECK(cell.amat_inv(6) == -3.0214007693_a);
      CHECK(cell.amat_inv(7) == 4.2423219287_a); 
      CHECK(cell.amat_inv(8) == -4.8560911922_a);

      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[0] == 0.121797_a);
      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[1] == -1.161093_a);
      CHECK(cell.crystal_to_cart(vec3d(0.2, -0.5, 0.867))[2] == -0.553419_a);

      CHECK(cell.cart_to_crystal(vec3d(0.66, -23.77, 2.72))[0] == -39.3396165136_a);
      CHECK(cell.cart_to_crystal(vec3d(0.66, -23.77, 2.72))[1] == 50.8091863243_a);
      CHECK(cell.cart_to_crystal(vec3d(0.66, -23.77, 2.72))[2] == -52.6483546581_a);

      CHECK(cell.contains(cell.crystal_to_cart(vec3d(0.5, 0.5, 0.5))));
      //This next one fails, this has to be checked.
      //CHECK(!cell.contains(cell.crystal_to_cart(vec3d(1.5, 0.5, 0.5))));
      CHECK(!cell.contains(cell.crystal_to_cart(vec3d(0.5, -0.1, 0.0))));
      CHECK(!cell.contains(cell.crystal_to_cart(vec3d(0.5, 0.5, -1.0))));
      
    }
  }
}
#endif
#endif
