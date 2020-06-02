/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__DIVERGENCE
#define OPERATIONS__DIVERGENCE

/*
 Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa, Alexey Karstev.

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

#include <basis/field.hpp>
#include <cassert>
namespace inq {
namespace operations {
	auto divergence(basis::field_set<basis::fourier_space, complex> const & ff){ // Divergence function for the field-set type defined in the Fourier space which return a filed 'diverg'
		basis::field<basis::fourier_space, complex> diverg(ff.basis());
		diverg = 0.0;
		for(int ix = 0; ix < ff.basis().sizes()[0]; ix++){		// Iterating over x-,y- and z- components of the input field 
			for(int iy = 0; iy < ff.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < ff.basis().sizes()[2]; iz++){
					auto gvec = ff.basis().gvector(ix, iy, iz);
					// Iterating over each vectorial components of input field-set and corresponding G-vector at (ix,iy,iz) point in the space
					for(int idir = 0; idir < 3 ; idir++) diverg.cubic()[ix][iy][iz] += complex(0.0, 1.0)*gvec[idir]*ff.cubic()[ix][iy][iz][idir]; 
				}
			}
		}
		return diverg;
	}
	auto divergence(basis::field_set<basis::real_space, complex> const & ff){	// Divergence function for the field-set type defined in the Real space which return a filed 'diverg_real'
		auto ff_fourier = operations::space::to_fourier(ff); 			// Tranform input field-set to Fourier space
		auto diverg_fourier = divergence(ff_fourier); 				// To calculate the divergence in Fourier space with the use of the above-defined function 'diverg'
		auto diverg_real = operations::space::to_real(diverg_fourier); 		// Transform output field to Real space
		return diverg_real;
		}
	auto divergence(basis::field_set<basis::real_space, double> const & ff){	// Divergence function for the field-set type defined in the Real space which return a filed 'diverg_real'
		auto ff_fourier = operations::space::to_fourier(ff.complex()); 		// Tranform input field-set to Fourier space
		auto diverg_fourier = divergence(ff_fourier); 				// To calculate the divergence in Fourier space with the use of the above-defined function 'diverg'
		auto diverg_real = operations::space::to_real(diverg_fourier); 		// Transform output field to Real space
		return diverg_real.real();						// Return a real part off the divergency in the real space
		}
}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vec3d.hpp>

        //Define test function 3
        auto f_analytic3 (inq::math::vec3d k , inq::math::vec3d r , int idir){
                using inq::math::vec3d;
		// Some random number to initialize the test function.
		//It ensures that the function will not clash with 0-values and will not give a false-positive test result
		auto f = inq::complex(30.30, 70.30);
		if ( idir == 0){f = 1.0 * exp(inq::complex(0.0, 1.0) * (k | r ));}
		if ( idir == 1){f = -2.3 * exp(inq::complex(0.0, 1.0) * (k | r ));}
		if ( idir == 2){f = 3.4 * exp(inq::complex(0.0, 1.0) * (k | r ));}
                return f;
        }

        //Define analytic form of the divergence of the test function 3
        auto d_analytic3 (inq::math::vec3d k , inq::math::vec3d r) {
                using inq::math::vec3d;
		// Some random number 
                auto g = inq::complex(70.30, 100.30);
                auto factor = inq::complex(0.0, 1.0)*exp(inq::complex(0.0,1.0)*(k | r ));
                g = factor * (1.0*k[0] - 2.3*k[1] + 3.4*k[2]);
                return g;
        }

        //Define test function 4
        auto f_analytic4 (inq::math::vec3d k , inq::math::vec3d r, int idir){
                using inq::math::vec3d;
		// Some random number 
                auto f = 30.90;
		if ( idir == 0){f = 1.0 * sin(k | r );}
		if ( idir == 1){f = -2.5 * cos(k | r );}
		if ( idir == 2){f = 3.3 * sin(k | r );}
                return f;
        }

        //Define analytic form of the divergence of the test function 4
        auto d_analytic4 (inq::math::vec3d k , inq::math::vec3d r) {
                using inq::math::vec3d;
		// Some random number 
                auto g =-80.03;
                g = (1.0 * k [0] * cos (k | r)) + (2.5 * k [1] * sin (k | r)) + (3.3 * k [2] * cos (k | r));
                return g;
        }


TEST_CASE("function operations::divergence", "[operations::divergence]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace operations;
	using math::vec3d;

	//UnitCell size
	double lx = 9;
	double ly = 12;
	double lz = 10;

	ions::geometry geo;
 	ions::UnitCell cell(vec3d(lx, 0.0, 0.0), vec3d(0.0, ly, 0.0), vec3d(0.0, 0.0, lz));

	basis::real_space rs(cell, input::basis::cutoff_energy(20.0));

	SECTION("Vectored plane-wave"){ 
		basis::field_set<basis::real_space, complex> f_test3(rs, 3);
	
		//Define k-vector for test function
		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) f_test3.cubic()[ix][iy][iz][idir] = f_analytic3 (kvec, vec, idir);
				}
			}
		}

		auto df_test3 = divergence(f_test3);

		double diff3 = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					diff3 += fabs(df_test3.cubic()[ix][iy][iz] - d_analytic3 (kvec, vec));
				}
			}
		}
		CHECK( diff3 < 1.0e-10 ); 
	}

	SECTION("Vectored real function"){

		basis::field_set<basis::real_space, double> f_test4(rs , 3);

		//Define k-vector for test function

		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) f_test4.cubic()[ix][iy][iz][idir] = f_analytic4 (kvec, vec, idir);
				}
			}
		}
		
		auto df_test4 = divergence(f_test4);

		double diff4 = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					diff4 += fabs(df_test4.cubic()[ix][iy][iz] - d_analytic4 (kvec, vec));
				}
			}
		}
		CHECK( diff4 < 1.0e-10 ); 
	}
}

#endif
#endif

