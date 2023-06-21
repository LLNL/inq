/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__LAPLACIAN
#define OPERATIONS__LAPLACIAN

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>

#include <gpu/run.hpp>
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <math/vector3.hpp>
#include <operations/space.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace operations {

template <typename FieldSetType, typename FactorType = double>
void laplacian_add(FieldSetType const & ff, FieldSetType & laplff, FactorType factor = 1.0, vector3<double, contravariant> const & gradcoeff = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;

	static_assert(std::is_same_v<typename FieldSetType::basis_type, basis::fourier_space>, "Only implemented for fourier_space");
			
	gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(), laplffcub = begin(laplff.hypercubic()), ffcub = begin(ff.hypercubic()), factor, gradcoeff]
					 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 auto lapl = factor*(-point_op.g2(ix, iy, iz) + dot(gradcoeff, point_op.gvector(ix, iy, iz)));
						 laplffcub[ix][iy][iz][ist] += lapl*ffcub[ix][iy][iz][ist];
					 });
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename FieldSetType, typename FactorType = double>
void laplacian_in_place(FieldSetType & ff, FactorType factor = 1.0, vector3<double, contravariant> const & gradcoeff = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;

	static_assert(std::is_same_v<typename FieldSetType::basis_type, basis::fourier_space>, "Only implemented for fourier_space");

	gpu::run(ff.set_part().local_size(), ff.basis().local_sizes()[2], ff.basis().local_sizes()[1], ff.basis().local_sizes()[0],
					 [point_op = ff.basis().point_op(),
						ffcub = begin(ff.hypercubic()), factor, gradcoeff] GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 auto lapl = factor*(-point_op.g2(ix, iy, iz) + dot(gradcoeff, point_op.gvector(ix, iy, iz)));						 
						 ffcub[ix][iy][iz][ist] = ffcub[ix][iy][iz][ist]*lapl;
					 });
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename FieldSetType, typename FactorType = double>
FieldSetType laplacian(FieldSetType const & ff, FactorType factor = 1.0, vector3<double, contravariant> const & gradcoeff = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;

	if constexpr(std::is_same_v<typename FieldSetType::basis_type, basis::real_space>) {
		return operations::space::to_real(operations::laplacian(operations::space::to_fourier(ff), factor));		
	} else {
		
		static_assert(std::is_same_v<typename FieldSetType::basis_type, basis::fourier_space>, "Only implemented for real or fourier_space");
		
		FieldSetType laplff(ff.skeleton());
		
		gpu::run(laplff.set_part().local_size(), laplff.basis().local_sizes()[2], laplff.basis().local_sizes()[1], laplff.basis().local_sizes()[0],
						 [point_op = ff.basis().point_op(), laplffcub = begin(laplff.hypercubic()), ffcub = begin(ff.hypercubic()), factor, gradcoeff]
						 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
						 auto lapl = factor*(-point_op.g2(ix, iy, iz) + dot(gradcoeff, point_op.gvector(ix, iy, iz)));
						 laplffcub[ix][iy][iz][ist] = lapl*ffcub[ix][iy][iz][ist];
						 });
		
		return laplff;
	}
}

}
}
#endif

#ifdef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST
#undef INQ_OPERATIONS_LAPLACIAN_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

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
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b);
		
		double factor = 0.673214;
	
		basis::real_space rs(box, /*spacing =*/ 0.49672941, basis_comm);

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
		systems::box box = systems::box::lattice({ll*1.0_b/sqrt(2), ll*1.0_b/sqrt(2), 0.0_b}, {-ll*1.0_b/sqrt(2), ll*1.0_b/sqrt(2), 0.0_b}, {0.0_b, 0.0_b, ll*1.0_b});
		
		double factor = 0.673214;
	
		basis::real_space rs(box, /*spacing =*/ 0.49672941, basis_comm);
		
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
		systems::box box = systems::box::lattice({0.0_b, ll*1.0_b, ll*1.0_b}, {ll*1.0_b, 0.0_b, ll*1.0_b}, {ll*1.0_b, ll*1.0_b, 0.0_b});
		
		double factor = 0.673214;
	
		basis::real_space rs(box, /*spacing =*/ 0.49672941, basis_comm);
		
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

