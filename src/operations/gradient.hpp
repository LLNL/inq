/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__GRADIENT
#define INQ__OPERATIONS__GRADIENT

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
#include <basis/fourier_space.hpp>
#include <operations/space.hpp>

#include <cassert>

namespace inq {
namespace operations {

	basis::field<basis::fourier_space, math::vector3<complex>> gradient(basis::field<basis::fourier_space, complex> const & ff){
		basis::field<basis::fourier_space, math::vector3<complex>> grad(ff.skeleton());

		gpu::run(grad.basis().local_sizes()[2], grad.basis().local_sizes()[1], grad.basis().local_sizes()[0],
						 [point_op = ff.basis().point_op(), gradcub = begin(grad.cubic()), ffcub = begin(ff.cubic())]
						 GPU_LAMBDA (auto iz, auto iy, auto ix){

							 auto gvec = point_op.gvector(ix, iy, iz);
							 for(int idir = 0; idir < 3; idir++) gradcub[ix][iy][iz][idir] = complex(0.0, 1.0)*gvec[idir]*ffcub[ix][iy][iz];
						 });

		return grad;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	basis::field_set<basis::fourier_space, math::vector3<complex>> gradient(basis::field_set<basis::fourier_space, complex> const & ff){
		basis::field_set<basis::fourier_space, math::vector3<complex>> grad(ff.skeleton());

		gpu::run(grad.set_part().local_size(), grad.basis().local_sizes()[2], grad.basis().local_sizes()[1], grad.basis().local_sizes()[0],
						 [point_op = ff.basis().point_op(), gradcub = begin(grad.cubic()), ffcub = begin(ff.cubic())]
						 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
							 
							 auto gvec = point_op.gvector(ix, iy, iz);
							 for(int idir = 0; idir < 3; idir++) gradcub[ix][iy][iz][ist][idir] = complex(0.0, 1.0)*gvec[idir]*ffcub[ix][iy][iz][ist];
						 });

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
		auto ff_fourier = operations::space::to_fourier(complex_field(ff));
		auto grad_fourier = gradient(ff_fourier);
		auto grad_real = operations::space::to_real(grad_fourier);
		return real_field(grad_real);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	auto gradient(basis::field_set<basis::real_space, complex> const & ff){
		auto ff_fourier = operations::space::to_fourier(ff);
		auto grad_fourier = gradient(ff_fourier);
		auto grad_real = operations::space::to_real(grad_fourier);
		return grad_real;
	}


	
}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_GRADIENT_UNIT_TEST
#undef INQ_OPERATIONS_GRADIENT_UNIT_TEST

#include <ions/geometry.hpp>

#include <catch2/catch.hpp>
#include <math/vector3.hpp>

auto f_analytic (inq::math::vector3<double> kk, inq::math::vector3<double> rr){
	return exp(inq::complex(0.0,1.0)*dot(kk, rr));
}

auto g_analytic (inq::math::vector3<double> kk , inq::math::vector3<double> rr) {
	inq::math::vector3<inq::complex> gg;
	auto factor = inq::complex(0.0, 1.0)*exp(inq::complex(0.0, 1.0)*dot(kk, rr));
	for(int idir = 0; idir < 3 ; idir++) gg[idir] = factor*kk[idir] ;
	return gg;
}

auto f_analytic2 (inq::math::vector3<double> kk, inq::math::vector3<double> rr){
	return sin(dot(kk, rr));
}

auto g_analytic2 (inq::math::vector3<double> kk , inq::math::vector3<double> rr) {
	inq::math::vector3<double> gg;
	for(int idir = 0; idir < 3 ; idir++) gg[idir] = kk[idir]*cos(dot(kk, rr));
	return gg;
}

TEST_CASE("function operations::gradient", "[operations::gradient]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace operations;
	using math::vector3;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

	//UnitCell size
	double lx = 9;
	double ly = 12;
	double lz = 10;
 	ions::UnitCell cell(vector3<double>(lx, 0.0, 0.0), vector3<double>(0.0, ly, 0.0), vector3<double>(0.0, 0.0, lz));

	SECTION("Plane-wave -- field"){ 

		basis::real_space rs(cell, input::basis::cutoff_energy(20.0), cart_comm);
	
		basis::field<basis::real_space, complex> f_test(rs);
	
		//Define k-vector for test function
		vector3<double> kvec = 2.0 * M_PI * vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					f_test.cubic()[ix][iy][iz] = f_analytic(kvec, vec);
				}
			}
		}

		auto g_test = gradient(f_test);

		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) diff += fabs(g_test.cubic()[ix][iy][iz][idir] - g_analytic (kvec, vec)[idir]);
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
				
		CHECK( diff < 1.0e-10 ); 
	}

	SECTION("Plane-wave -- field_set"){

		basis::real_space rs(cell, input::basis::cutoff_energy(20.0), basis_comm);
		
		basis::field_set<basis::real_space, complex> f_test(rs, 13, cart_comm);
	
		//Define k-vector for test function
		vector3<double> kvec = 2.0*M_PI*vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int ist = 0; ist < f_test.local_set_size(); ist++){					
						f_test.cubic()[ix][iy][iz][ist] = double(ist)*f_analytic(kvec, vec);
					}
				}
			}
		}

		auto g_test = gradient(f_test);

		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int ist = 0; ist < f_test.local_set_size(); ist++){
						for(int idir = 0; idir < 3 ; idir++) diff += fabs(g_test.cubic()[ix][iy][iz][ist][idir] - double(ist)*g_analytic(kvec, vec)[idir]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		CHECK( diff < 1.0e-8 ); 
	}
		
	SECTION("Real function"){

		basis::real_space rs(cell, input::basis::cutoff_energy(20.0), cart_comm);
		
		basis::field<basis::real_space, double> f_test2(rs);
	
		vector3<double> kvec = 2.0 * M_PI * vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					f_test2.cubic()[ix][iy][iz] = f_analytic2(kvec, vec);
				}
			}
		}

		auto g_test2 = gradient(f_test2);
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					for(int idir = 0; idir < 3 ; idir++) diff += fabs(g_test2.cubic()[ix][iy][iz][idir] - g_analytic2(kvec, vec)[idir]);
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
				
		CHECK( diff < 1.0e-10 ); 
	}
}

#endif
#endif

