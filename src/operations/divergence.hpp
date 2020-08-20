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

	auto point_op = ff.basis().point_op();
	
	for(int ix = 0; ix < ff.basis().sizes()[0]; ix++){
		for(int iy = 0; iy < ff.basis().sizes()[1]; iy++){
			for(int iz = 0; iz < ff.basis().sizes()[2]; iz++){

				auto gvec = point_op.gvector(ix, iy, iz);
				complex div = 0.0;

				for(int idir = 0; idir < 3 ; idir++) div += complex(0.0, 1.0)*gvec[idir]*ff.cubic()[ix][iy][iz][idir]; 
				diverg.cubic()[ix][iy][iz] = div;
			}
		}
	}
	
	return diverg;
}

auto divergence(basis::field_set<basis::real_space, complex> const & ff){
	auto ff_fourier = operations::space::to_fourier(ff); 			
	auto diverg_fourier = divergence(ff_fourier); 				
	auto diverg_real = operations::space::to_real(diverg_fourier);
	return diverg_real;
}

auto divergence(basis::field_set<basis::real_space, double> const & ff){
	auto ff_fourier = operations::space::to_fourier(ff.complex());	
	auto diverg_fourier = divergence(ff_fourier); 
	auto diverg_real = operations::space::to_real(diverg_fourier);
	return diverg_real.real();
}

}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vector3.hpp>

		//Define test function 3
		auto vectorial_complex_plane_wave (inq::math::vec3d k , inq::math::vec3d r ){
				using inq::math::vec3d;
				std::array<inq::complex, 3> f;
				f[0] = 1.0*exp(inq::complex(0.0, 1.0)*(k | r ));
				f[1] = -2.3*exp(inq::complex(0.0, 1.0)*(k | r ));
				f[2] = 3.4*exp(inq::complex(0.0, 1.0)*(k | r ));
				return f;
				}

		//Define analytic form of the divergence of the test vectorial_complex_plane_wave
		auto d_vectorial_complex_plane_wave (inq::math::vec3d k , inq::math::vec3d r) {
				using inq::math::vec3d;
				// Some random number 
				auto factor = inq::complex(0.0, 1.0)*exp(inq::complex(0.0,1.0)*(k | r ));
				return factor*(1.0*k[0] - 2.3*k[1] + 3.4*k[2]);
		}

		//Define test function 4
		auto vectorial_real_wave (inq::math::vec3d k , inq::math::vec3d r){
				using inq::math::vec3d;
				std::array<double, 3> f;
				// Some random number 
				f[0] = 1.0*sin(k | r );
				f[1] = -2.5*cos(k | r );
				f[2] = 3.3*sin(k | r );
				return f;
		}

		//Define analytic form of the divergence of the test function 4
		auto d_vectorial_real_wave (inq::math::vec3d k , inq::math::vec3d r) {
				using inq::math::vec3d;
				return (1.0*k[0]*cos(k | r)) + (2.5*k[1]*sin(k | r)) + (3.3*k[2]*cos(k | r));
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
		basis::field_set<basis::real_space, complex> vectorial_complex_field(rs, 3);
	
		//Define k-vector for test function
		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) vectorial_complex_field.cubic()[ix][iy][iz][idir] = vectorial_complex_plane_wave (kvec, vec)[idir];
				}
			}
		}

		auto d_vectorial_complex_field = divergence(vectorial_complex_field);

		double diff3 = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					diff3 += fabs(d_vectorial_complex_field.cubic()[ix][iy][iz] - d_vectorial_complex_plane_wave (kvec, vec));
				}
			}
		}
		CHECK( diff3 < 1.0e-10 ); 
	}

	SECTION("Vectored real function"){

		basis::field_set<basis::real_space, double> vectorial_real_field(rs , 3);

		//Define k-vector for test function

		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) vectorial_real_field.cubic()[ix][iy][iz][idir] = vectorial_real_wave (kvec, vec)[idir];
				}
			}
		}
		
		auto d_vectorial_real_field = divergence(vectorial_real_field);

		double diff4 = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					diff4 += fabs(d_vectorial_real_field.cubic()[ix][iy][iz] - d_vectorial_real_wave (kvec, vec));
				}
			}
		}
		CHECK( diff4 < 1.0e-10 ); 
	}
}

#endif
#endif

