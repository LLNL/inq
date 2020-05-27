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

	auto gradient(basis::field<basis::fourier_space, complex> const & ff){
		basis::field_set<basis::fourier_space, complex> grad(ff.basis(), 3);

		for(int ix = 0; ix < ff.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < ff.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < ff.basis().sizes()[2]; iz++){
					auto gvec = ff.basis().gvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) grad.cubic()[ix][iy][iz][idir] = complex(0.0, 1.0)*gvec[idir]*ff.cubic()[ix][iy][iz]; 
				}
			}
		}
		return grad;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	auto gradient(basis::field<basis::real_space, complex> const & ff){
		auto ff_fourier = operations::space::to_fourier(ff);
		auto grad_fourier = gradient(ff_fourier);
		auto grad_real = operations::space::to_real(grad_fourier);
		return grad_real;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	auto gradient(basis::field<basis::real_space, double> const & ff){
		auto ff_fourier = operations::space::to_fourier(ff);
		auto grad_fourier = gradient(ff_fourier);
		auto grad_real = operations::space::to_real(grad_fourier);
		return grad_real.real();
	}
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vec3d.hpp>
#include <math/vector3.hpp>

auto f_analytic (math::vec3d kk, math::vec3d rr){
	return exp(complex(0.0,1.0)*(kk|rr));
}

auto g_analytic (math::vec3d kk , math::vec3d rr) {
	math::vector3<complex> gg;
	complex factor = complex(0.0, 1.0)*exp(complex(0.0, 1.0)*(kk|rr));
	for(int idir = 0; idir < 3 ; idir++) gg[idir] = factor*kk[idir] ;
	return gg;
}

auto f_analytic2 (math::vec3d kk, math::vec3d rr){
	return sin(kk|rr);
}

auto g_analytic2 (math::vec3d kk , math::vec3d rr) {
	math::vec3d gg;
	for(int idir = 0; idir < 3 ; idir++) gg[idir] = kk[idir]*cos(kk|rr);
	return gg;
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

		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					f_test.cubic()[ix][iy][iz] = f_analytic(kvec, vec);
				}
			}
		}

		auto g_test = gradient(f_test);

		double diff = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) diff += fabs(g_test.cubic()[ix][iy][iz][idir] - g_analytic (kvec, vec)[idir]);
				}
			}
		}
		
		CHECK( diff < 1.0e-10 ); 
	}
	
	SECTION("Real function"){

		basis::field<basis::real_space, double> f_test2(rs);
	
		vec3d kvec = 2.0 * M_PI * vec3d(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					f_test2.cubic()[ix][iy][iz] = f_analytic2(kvec, vec);
				}
			}
		}

		auto g_test2 = gradient(f_test2);
		double diff2 = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) diff2 += fabs(g_test2.cubic()[ix][iy][iz][idir] - g_analytic2 (kvec, vec)[idir]);
				}
			}
		}
		
		CHECK( diff2 < 1.0e-10 ); 
	}
}

#endif
#endif

