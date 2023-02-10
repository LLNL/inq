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

#include <cassert>

#include <gpu/run.hpp>
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <math/vector3.hpp>
#include <operations/space.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace operations {


template <template<typename, typename> class SetType, typename FactorType = double>
void laplacian_add(SetType<basis::fourier_space, complex> const & ff, SetType<basis::fourier_space, complex>& laplff, FactorType factor = 1.0, vector3<double, contravariant> const & gradcoeff = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;
		
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(), laplffcub = begin(laplff.hypercubic()), ffcub = begin(ff.hypercubic()), factor, gradcoeff]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 auto lapl = factor*(-point_op.g2(ix, iy, iz) + dot(gradcoeff, point_op.gvector(ix, iy, iz)));
						 laplffcub[ix][iy][iz][ist] += lapl*ffcub[ix][iy][iz][ist];
					 });
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <template<typename, typename> class SetType, typename FactorType = double>
void laplacian_in_place(SetType<basis::fourier_space, complex>& ff, FactorType factor = 1.0, vector3<double, contravariant> const & gradcoeff = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;
		
	gpu::run(ff.set_part().local_size(), ff.basis().local_sizes()[2], ff.basis().local_sizes()[1], ff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						ffcub = begin(ff.hypercubic()), factor, gradcoeff] GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 auto lapl = factor*(-point_op.g2(ix, iy, iz) + dot(gradcoeff, point_op.gvector(ix, iy, iz)));						 
						 ffcub[ix][iy][iz][ist] = ffcub[ix][iy][iz][ist]*lapl;
					 });
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <template<typename, typename> class SetType, typename FactorType = double>
SetType<basis::fourier_space, complex> laplacian(SetType<basis::fourier_space, complex> const & ff, FactorType factor = 1.0, vector3<double, contravariant> const & gradcoeff = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;
	
	SetType<basis::fourier_space, complex> laplff(ff.skeleton());
	
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(), laplffcub = begin(laplff.hypercubic()), ffcub = begin(ff.hypercubic()), factor, gradcoeff]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 auto lapl = factor*(-point_op.g2(ix, iy, iz) + dot(gradcoeff, point_op.gvector(ix, iy, iz)));
						 laplffcub[ix][iy][iz][ist] = lapl*ffcub[ix][iy][iz][ist];
					 });

	return laplff;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <template<typename, typename> class SetType, typename FactorType = double>
auto laplacian(SetType<basis::real_space, complex> const & ff, FactorType factor = 1.0){

	return operations::space::to_real(operations::laplacian(operations::space::to_fourier(ff), factor));
}

}
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST
#undef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("function operations::laplacian", "[operations::laplacian]") {

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using namespace operations;
	using Catch::Approx;
	
	parallel::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	auto set_comm = basis::set_subcomm(cart_comm);
	auto basis_comm = basis::basis_subcomm(cart_comm);	

	SECTION("Plane-wave -- field_set"){
		
		double lx = 9;
		double ly = 12;
		double lz = 10;
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);
		
		double factor = 0.673214;
	
		basis::real_space rs(box, basis_comm);
		
		basis::field_set<basis::real_space, complex> func(rs, 13, cart_comm);
	
		//Define k-vector for test function
		auto kvec = 2.0*M_PI*vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);
		
		auto ff = [] (auto & kk, auto & rr){
			return exp(inq::complex(0.0, 1.0)*dot(kk, rr));
		};
		
		auto laplff = [ff] (auto & kk, auto & rr) {
			return -norm(kk)*ff(kk, rr);
		};
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.hypercubic()[ix][iy][iz][ist] = (ist + 1.0)*ff(kvec, vec);
				}
			}
		}

		auto lapl = operations::laplacian(func, factor);
		auto func_fs = operations::space::to_fourier(func);
		operations::laplacian_in_place(func_fs, factor);
		auto lapl_in_place = operations::space::to_real(func_fs);
		operations::laplacian_add(operations::space::to_fourier(func), func_fs, factor);
		auto lapl_add = operations::space::to_real(func_fs);
		
		double diff = 0.0;
		double diff_in_place = 0.0;
		double diff_add = 0.0;		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(lapl.hypercubic()[ix][iy][iz][ist] - anvalue);
						diff_in_place += fabs(lapl_in_place.hypercubic()[ix][iy][iz][ist] - anvalue);
						diff_add += fabs(lapl_add.hypercubic()[ix][iy][iz][ist] - 2.0*anvalue);						
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_in_place, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_add, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-8) ;
		CHECK(diff_in_place < 1.0e-8);
		CHECK(diff_add < 5.0e-7);
		
	}
	
	SECTION("Plane-wave -- rotated"){
		
		double ll = 9;
		systems::box box = systems::box::lattice({ll*1.0_b/sqrt(2), ll*1.0_b/sqrt(2), 0.0_b}, {-ll*1.0_b/sqrt(2), ll*1.0_b/sqrt(2), 0.0_b}, {0.0_b, 0.0_b, ll*1.0_b}).cutoff_energy(20.0_Ha);
		
		double factor = 0.673214;
	
		basis::real_space rs(box, basis_comm);

		CHECK(rs.cell().volume() == ll*ll*ll);
		
		basis::field_set<basis::real_space, complex> func(rs, 13, cart_comm);
	
		auto kvec = 2.0*M_PI*vector3<double>(sqrt(2.0)/ll, sqrt(2.0)/ll, 0.0);
		
		auto ff = [] (auto & kk, auto & rr){
			return exp(inq::complex(0.0, 1.0)*dot(kk, rr));
		};
		
		auto laplff = [ff] (auto & kk, auto & rr) {
			return -norm(kk)*ff(kk, rr);
		};

		CHECK(norm(kvec) == Approx(rs.cell().metric().norm(rs.cell().metric().to_covariant(kvec))));
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.hypercubic()[ix][iy][iz][ist] = (ist + 1.0)*ff(kvec, vec);
				}
			}
		}

		auto lapl = operations::laplacian(func, factor);
		auto func_fs = operations::space::to_fourier(func);
		operations::laplacian_in_place(func_fs, factor);
		auto lapl_in_place = operations::space::to_real(func_fs);
		operations::laplacian_add(operations::space::to_fourier(func), func_fs, factor);
		auto lapl_add = operations::space::to_real(func_fs);
		
		double diff = 0.0;
		double diff_in_place = 0.0;
		double diff_add = 0.0;		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(lapl.hypercubic()[ix][iy][iz][ist] - anvalue);
						diff_in_place += fabs(lapl_in_place.hypercubic()[ix][iy][iz][ist] - anvalue);
						diff_add += fabs(lapl_add.hypercubic()[ix][iy][iz][ist] - 2.0*anvalue);						
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_in_place, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_add, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-8) ;
		CHECK(diff_in_place < 1.0e-8);
		CHECK(diff_add < 5.0e-7);
		
	}
	
	SECTION("Plane-wave -- non-orthogonal"){
		
		double ll = 5.89;
		systems::box box = systems::box::lattice({0.0_b, ll*1.0_b, ll*1.0_b}, {ll*1.0_b, 0.0_b, ll*1.0_b}, {ll*1.0_b, ll*1.0_b, 0.0_b}).cutoff_energy(20.0_Ha);
		
		double factor = 0.673214;
	
		basis::real_space rs(box, basis_comm);

		basis::field_set<basis::real_space, complex> func(rs, 13, cart_comm);
	
		auto kvec = rs.cell().metric().to_cartesian(2.0*rs.cell().reciprocal(2) + 3.0*rs.cell().reciprocal(2) - 1.0*rs.cell().reciprocal(2));
		
		auto ff = [] (auto & kk, auto & rr){
			return exp(inq::complex(0.0, 1.0)*dot(kk, rr));
		};
		
		auto laplff = [ff] (auto & kk, auto & rr) {
			return -norm(kk)*ff(kk, rr);
		};

		CHECK(norm(kvec) == Approx(rs.cell().metric().norm(rs.cell().metric().to_covariant(kvec))));
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.hypercubic()[ix][iy][iz][ist] = (ist + 1.0)*ff(kvec, vec);
				}
			}
		}

		auto lapl = operations::laplacian(func, factor);
		auto func_fs = operations::space::to_fourier(func);
		operations::laplacian_in_place(func_fs, factor);
		auto lapl_in_place = operations::space::to_real(func_fs);
		operations::laplacian_add(operations::space::to_fourier(func), func_fs, factor);
		auto lapl_add = operations::space::to_real(func_fs);
		
		double diff = 0.0;
		double diff_in_place = 0.0;
		double diff_add = 0.0;		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(lapl.hypercubic()[ix][iy][iz][ist] - anvalue);
						diff_in_place += fabs(lapl_in_place.hypercubic()[ix][iy][iz][ist] - anvalue);
						diff_add += fabs(lapl_add.hypercubic()[ix][iy][iz][ist] - 2.0*anvalue);						
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_in_place, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&diff_add, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-7) ;
		CHECK(diff_in_place < 1.0e-7);
		CHECK(diff_add < 5.0e-7);
		
	}
	
}
#endif

