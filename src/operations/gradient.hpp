/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__GRADIENT
#define OPERATIONS__GRADIENT

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

namespace operations {
	auto gradient(basis::field<basis::fourier_space, complex> const & ff){		// Gradient function for the field type defined in the Fourier space which return a filed-set 'grad'
		basis::field_set<basis::fourier_space, complex> grad(ff.basis(), 3);
		for(int ix = 0; ix < ff.basis().sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < ff.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < ff.basis().sizes()[2]; iz++){
					auto gvec = ff.basis().gvector(ix, iy, iz);
					// Iterating over each vectorial components of grad field-set and corresponding G-vector at (ix,iy,iz) point in the space
					for(int idir = 0; idir < 3 ; idir++) grad.cubic()[ix][iy][iz][idir] = complex(0.0, 1.0)*gvec[idir]*ff.cubic()[ix][iy][iz]; 
				}
			}
		}
		return grad;
	}
	auto gradient(basis::field<basis::real_space, complex> const & ff){		// Gradient function for the field type (complex) defined in the Real space which return a filed-set 'grad_real'
		auto ff_fourier = operations::space::to_fourier(ff); 			// Tranform input field to Fourier space
		auto grad_fourier = gradient(ff_fourier); 				// To calculate the gradient in Fourier space with the use of the above-defined function 'gradient'
		auto grad_real = operations::space::to_real(grad_fourier); 		// Transform output field-set to Real space (complex)
		return grad_real;
	}
	auto gradient(basis::field<basis::real_space, double> const & ff){		// Gradient function for the field type (double) defined in the Real space which return a filed-set 'grad_real'
		auto ff_fourier = operations::space::to_fourier(ff); 			// Tranform input field to Fourier space
		auto grad_fourier = gradient(ff_fourier); 				// To calculate the gradient in Fourier space with the use of the above-defined function 'gradient'
		auto grad_real = operations::space::to_real(grad_fourier); 		// Transform output field-set to Real space (double)
		return grad_real;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vec3d.hpp>

	//Define test function 1
	complex f_analytic (math::vec3d k , math::vec3d r){
		using math::vec3d;
		complex f;
		f = exp(complex(0.0,1.0)*(k | r ));
		return f;
	}

	//Define analytic form of the gradient of the test function 1
	complex* g_analytic (math::vec3d k , math::vec3d r) {
		using math::vec3d;
		static complex g[3];
		complex factor = complex(0.0, 1.0)*exp(complex(0.0,1.0)*(k | r ));
		for(int idir = 0; idir < 3 ; idir++) g [idir] = factor * k [idir] ;
		return g;
	}

	//Define test function 2
	double f_analytic2 (math::vec3d k , math::vec3d r){
		using math::vec3d;
		double f;
		f = sin(k | r );
		return f;
	}

	//Define analytic form of the gradient of the test function 1
	double* g_analytic2 (math::vec3d k , math::vec3d r) {
		using math::vec3d;
		static double g[3];
		for(int idir = 0; idir < 3 ; idir++) g [idir] = k [idir] * cos (k | r);
		return g;
	}

TEST_CASE("function operations::gradient", "[operations::gradient]") {

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

	SECTION("Plane-wave"){ 
		basis::field<basis::real_space, complex> f_test(rs);
	
		//Define k-vector for test function
		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					f_test.cubic()[ix][iy][iz] = f_analytic (kvec, vec);
				}
			}
		}

		basis::field_set<basis::real_space, complex> g_test(rs, 3);
		g_test = gradient(f_test);

		double diff = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) diff += abs(g_test.cubic()[ix][iy][iz][idir] - g_analytic (kvec, vec)[idir]);
				}
			}
		}
		CHECK( diff < 1.0e-10 ); 
	}
	
	SECTION("Real function"){

		basis::field<basis::real_space, double> f_test2(rs);
	
		//Define k-vector for test function
		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					f_test2.cubic()[ix][iy][iz] = f_analytic2 (kvec, vec);
				}
			}
		}

		basis::field_set<basis::real_space, complex> g_test2(rs, 3);
		g_test2 = gradient(f_test2);

		double diff2 = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) diff2 += abs(g_test2.cubic()[ix][iy][iz][idir] - g_analytic2 (kvec, vec)[idir]);
				}
			}
		}
		CHECK( diff2 < 1.0e-10 ); 
	}
}

#endif
#endif

