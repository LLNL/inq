/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__LAPLACIAN
#define OPERATIONS__LAPLACIAN

/*
 Copyright (C) 2020-2021 Xavier Andrade, Alfredo A. Correa.

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

#include <gpu/run.hpp>
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <operations/space.hpp>

#include <multi/adaptors/fftw.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/cufft.hpp>
#endif

#include <cassert>

#include <utils/profiling.hpp>

namespace inq {
namespace operations {
	
void laplacian_add(basis::field_set<basis::fourier_space, complex> const & ff, basis::field_set<basis::fourier_space, complex>& laplff){

	CALI_CXX_MARK_FUNCTION;
		
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						laplffcub = begin(laplff.cubic()),
						ffcub = begin(ff.cubic())]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 double lapl = -0.5*(-point_op.g2(ix, iy, iz));
						 laplffcub[ix][iy][iz][ist] += lapl*ffcub[ix][iy][iz][ist];
					 });
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void laplacian_in_place(basis::field_set<basis::fourier_space, complex>& ff){

	CALI_CXX_MARK_FUNCTION;
		
	gpu::run(ff.set_part().local_size(), ff.basis().local_sizes()[2], ff.basis().local_sizes()[1], ff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						ffcub = begin(ff.cubic())] GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 double lapl = -0.5*(-point_op.g2(ix, iy, iz));
						 ffcub[ix][iy][iz][ist] = ffcub[ix][iy][iz][ist]*lapl;
					 });
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

basis::field_set<basis::fourier_space, complex> laplacian(basis::field_set<basis::fourier_space, complex> const & ff){

	CALI_CXX_MARK_FUNCTION;
	
	basis::field_set<basis::fourier_space, complex> laplff(ff.skeleton());
	
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						laplffcub = begin(laplff.cubic()),
						ffcub = begin(ff.cubic())]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 double lapl = -0.5*(-point_op.g2(ix, iy, iz));
						 laplffcub[ix][iy][iz][ist] = lapl*ffcub[ix][iy][iz][ist];
					 });

	return laplff;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto laplacian(basis::field_set<basis::real_space, complex> const & ff){

	return operations::space::to_real(operations::laplacian(operations::space::to_fourier(ff)));
}

}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST
#undef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST

#include <catch2/catch.hpp>

auto ff(inq::math::vector3<double> const & kk, inq::math::vector3<double> const & rr){
	return exp(inq::complex(0.0,1.0)*dot(kk, rr));
}

auto laplff(inq::math::vector3<double> const & kk, inq::math::vector3<double> const & rr) {
	return -dot(kk, kk)*ff(kk, rr);
}

TEST_CASE("function operations::gradient", "[operations::gradient]") {

	using namespace inq;
	using namespace inq::magnitude;	
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

	double factor = -0.5;
	
	SECTION("Plane-wave -- field_set"){

		basis::real_space rs(cell, input::basis::cutoff_energy(20.0_Ha), basis_comm);
		
		basis::field_set<basis::real_space, complex> func(rs, 13, cart_comm);
	
		//Define k-vector for test function
		vector3<double> kvec = 2.0*M_PI*vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.cubic()[ix][iy][iz][ist] = (ist + 1.0)*ff(kvec, vec);
				}
			}
		}

		auto lapl = operations::laplacian(func);
		auto func_fs = operations::space::to_fourier(func);
		operations::laplacian_in_place(func_fs);
		auto lapl_in_place = operations::space::to_real(func_fs);
		
		double diff = 0.0;
		double diff_in_place = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(lapl.cubic()[ix][iy][iz][ist] - anvalue);
						diff_in_place += fabs(lapl.cubic()[ix][iy][iz][ist] - anvalue);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_in_place, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-8) ;
		CHECK(diff_in_place < 1.0e-8);
				
	}

}

#endif

#endif

