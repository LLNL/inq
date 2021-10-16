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
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <operations/space.hpp>

#include <cassert>
namespace inq {
namespace operations {

basis::field<basis::fourier_space, complex> divergence(basis::field<basis::fourier_space, math::vector3<complex>> const & ff){ // Divergence function for the field-set type defined in the Fourier space which return a filed 'diverg'
	
	basis::field<basis::fourier_space, complex> diverg(ff.basis());

	gpu::run(diverg.basis().local_sizes()[2], diverg.basis().local_sizes()[1], diverg.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(), divergcub = begin(diverg.cubic()), ffcub = begin(ff.cubic())]
					 GPU_LAMBDA (auto iz, auto iy, auto ix){

						 auto gvec = point_op.gvector(ix, iy, iz);
						 complex div = 0.0;
						 for(int idir = 0; idir < 3 ; idir++) div += complex(0.0, 1.0)*gvec[idir]*ffcub[ix][iy][iz][idir]; 
						 divergcub[ix][iy][iz] = div;
					 });
	
	return diverg;
}

auto divergence(basis::field<basis::real_space, math::vector3<complex>> const & ff){
	auto ff_fourier = operations::space::to_fourier(ff); 			
	auto diverg_fourier = divergence(ff_fourier); 				
	auto diverg_real = operations::space::to_real(diverg_fourier);
	return diverg_real;
}

auto divergence(basis::field<basis::real_space, math::vector3<double>> const & ff){
	auto ff_fourier = operations::space::to_fourier(complex_field(ff));
	auto diverg_fourier = divergence(ff_fourier); 
	auto diverg_real = operations::space::to_real(diverg_fourier);
	return real_field(diverg_real);
}

}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_DIVERGENCE_UNIT_TEST
#undef INQ_OPERATIONS_DIVERGENCE_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/vector3.hpp>

auto vectorial_complex_plane_wave (inq::math::vector3<double> k , inq::math::vector3<double> r ){
	std::array<inq::complex, 3> f;
	f[0] = 1.0*exp(inq::complex(0.0, 1.0)*dot(k, r));
	f[1] = -2.3*exp(inq::complex(0.0, 1.0)*dot(k, r));
	f[2] = 3.4*exp(inq::complex(0.0, 1.0)*dot(k, r));
	return f;
}

auto d_vectorial_complex_plane_wave (inq::math::vector3<double> k , inq::math::vector3<double> r) {
	auto factor = inq::complex(0.0, 1.0)*exp(inq::complex(0.0,1.0)*dot(k, r));
	return factor*(1.0*k[0] - 2.3*k[1] + 3.4*k[2]);
}

auto vectorial_real_wave (inq::math::vector3<double> k , inq::math::vector3<double> r){
	std::array<double, 3> f;

	f[0] = 1.0*sin(dot(k, r));
	f[1] = -2.5*cos(dot(k, r));
	f[2] = 3.3*sin(dot(k,r));
	return f;
}

auto d_vectorial_real_wave (inq::math::vector3<double> k , inq::math::vector3<double> r) {
	return 1.0*k[0]*cos(dot(k, r)) + 2.5*k[1]*sin(dot(k, r)) + 3.3*k[2]*cos(dot(k, r));
}

TEST_CASE("function operations::divergence", "[operations::divergence]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using namespace operations;
	using math::vector3;

	//UnitCell size
	double lx = 9;
	double ly = 12;
	double lz = 10;
	systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	//	auto set_comm = cart_comm.axis(0);
	//	auto basis_comm = cart_comm.axis(1);

	basis::real_space rs(box, cart_comm);

	SECTION("Vectored plane-wave"){ 
		basis::field<basis::real_space, math::vector3<complex>> vectorial_complex_field(rs);
	
		//Define k-vector for test function
		vector3<double> kvec = 2.0 * M_PI * vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) vectorial_complex_field.cubic()[ix][iy][iz][idir] = vectorial_complex_plane_wave (kvec, vec)[idir];
				}
			}
		}

		auto d_vectorial_complex_field = divergence(vectorial_complex_field);

		double diff3 = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector(ix, iy, iz);
					diff3 += fabs(d_vectorial_complex_field.cubic()[ix][iy][iz] - d_vectorial_complex_plane_wave (kvec, vec));
				}
			}
		}
		CHECK( diff3 < 1.0e-10 ); 
	}

	SECTION("Vectored real function"){

		basis::field<basis::real_space, math::vector3<double>> vectorial_real_field(rs);

		//Define k-vector for test function

		vector3<double> kvec = 2.0 * M_PI * vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field 
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) vectorial_real_field.cubic()[ix][iy][iz][idir] = vectorial_real_wave (kvec, vec)[idir];
				}
			}
		}
		
		auto d_vectorial_real_field = divergence(vectorial_real_field);

		double diff4 = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){ 			// Iterating over each x-,y- and z- components of the input field-set 
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector(ix, iy, iz);
					diff4 += fabs(d_vectorial_real_field.cubic()[ix][iy][iz] - d_vectorial_real_wave (kvec, vec));
				}
			}
		}
		CHECK( diff4 < 1.0e-10 ); 
	}
}

#endif
#endif

