/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__DIVERGENCE
#define OPERATIONS__DIVERGENCE

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

template <typename FieldSetType, typename FactorType = double,
					typename ResultType = typename FieldSetType::template template_type<typename FieldSetType::basis_type, typename FieldSetType::element_type::element_type>>
ResultType divergence(FieldSetType const & ff, double factor = 1.0, vector3<double, covariant> const & shift = {0.0, 0.0, 0.0}){

	CALI_CXX_MARK_FUNCTION;

	if constexpr(std::is_same_v<typename FieldSetType::basis_type, basis::real_space>) {
		if constexpr(std::is_same_v<typename FieldSetType::element_type::element_type, double>) {		
			return real_field(operations::space::to_real(operations::divergence(operations::space::to_fourier(complex_field(ff)), factor, shift)));
		} else {
			return operations::space::to_real(operations::divergence(operations::space::to_fourier(ff), factor, shift));
		}
	} else {
		
		static_assert(std::is_same_v<typename FieldSetType::basis_type, basis::fourier_space>, "Only implemented for real or fourier_space");
		ResultType divff(ff.skeleton());

		gpu::run(divff.set_part().local_size(), divff.basis().local_sizes()[2], divff.basis().local_sizes()[1], divff.basis().local_sizes()[0],
						 [point_op = ff.basis().point_op(), divffcub = begin(divff.hypercubic()), ffcub = begin(ff.hypercubic()), factor, shift]
						 GPU_LAMBDA (auto ist, auto iz, auto iy, auto ix){
							 auto gvec = point_op.gvector(ix, iy, iz) + shift; 
							 divffcub[ix][iy][iz][ist] = factor*complex(0.0, 1.0)*point_op.metric().dot(gvec, ffcub[ix][iy][iz][ist]);
						 });
		
		return divff;
	}
}

}
}
#endif

#ifdef INQ_OPERATIONS_DIVERGENCE_UNIT_TEST
#undef INQ_OPERATIONS_DIVERGENCE_UNIT_TEST

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
		
		basis::field_set<basis::real_space, vector3<complex, cartesian>> func(rs, 13, cart_comm);
	
		//Define k-vector for test function
		auto kvec = 2.0*M_PI*vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);
		
		auto ff = [] (auto & kk, auto & rr){
			return exp(inq::complex(0.0, 1.0)*dot(kk, rr));
		};
		
		auto gradff = [ff] (auto & kk, auto & rr) {
			return complex(0.0, 1.0)*kk*ff(kk, rr);
		};

		auto laplff = [ff] (auto & kk, auto & rr) {
			return -norm(kk)*ff(kk, rr);
		};
				
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.hypercubic()[ix][iy][iz][ist] = (ist + 1.0)*gradff(kvec, vec);
				}
			}
		}

		auto div = operations::divergence(func, factor);
		
		auto diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(div.hypercubic()[ix][iy][iz][ist] - anvalue);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-8) ;
		
	}
	
	SECTION("Plane-wave -- rotated"){
		
		double ll = 9;
		systems::box box = systems::box::lattice({ll*1.0_b/sqrt(2), ll*1.0_b/sqrt(2), 0.0_b}, {-ll*1.0_b/sqrt(2), ll*1.0_b/sqrt(2), 0.0_b}, {0.0_b, 0.0_b, ll*1.0_b});
		
		double factor = 0.673214;
	
		basis::real_space rs(box, /*spacing =*/ 0.49672941, basis_comm);

		CHECK(rs.cell().volume() == ll*ll*ll);
		
		basis::field_set<basis::real_space, vector3<complex, cartesian>> func(rs, 13, cart_comm);
	
		auto kvec = 2.0*M_PI*vector3<double>(sqrt(2.0)/ll, sqrt(2.0)/ll, 0.0);
		
		auto ff = [] (auto & kk, auto & rr){
			return exp(inq::complex(0.0, 1.0)*dot(kk, rr));
		};
		
		auto gradff = [ff] (auto & kk, auto & rr) {
			return complex(0.0, 1.0)*kk*ff(kk, rr);			
		};

		auto laplff = [ff] (auto & kk, auto & rr) {
			return -norm(kk)*ff(kk, rr);
		};
		
		CHECK(norm(kvec) == Approx(rs.cell().metric().norm(rs.cell().metric().to_covariant(kvec))));
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.hypercubic()[ix][iy][iz][ist] = (ist + 1.0)*gradff(kvec, vec);
				}
			}
		}

		auto div = operations::divergence(func, factor);
		
		auto diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(div.hypercubic()[ix][iy][iz][ist] - anvalue);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-8);
	}
	
	SECTION("Plane-wave -- non-orthogonal"){
		
		double ll = 5.89;
		systems::box box = systems::box::lattice({0.0_b, ll*1.0_b, ll*1.0_b}, {ll*1.0_b, 0.0_b, ll*1.0_b}, {ll*1.0_b, ll*1.0_b, 0.0_b});
		
		double factor = 0.673214;

		basis::real_space rs(box, /*spacing =*/ 0.49672941, basis_comm);

		basis::field_set<basis::real_space, vector3<complex, cartesian>> func(rs, 13, cart_comm);
	
		auto kvec = rs.cell().metric().to_cartesian(2.0*rs.cell().reciprocal(2) + 3.0*rs.cell().reciprocal(2) - 1.0*rs.cell().reciprocal(2));
		
		auto ff = [] (auto & kk, auto & rr){
			return exp(inq::complex(0.0, 1.0)*dot(kk, rr));
		};
		
		auto gradff = [ff] (auto & kk, auto & rr) {
			return complex(0.0, 1.0)*kk*ff(kk, rr);
		};

		auto laplff = [ff] (auto & kk, auto & rr) {
			return -norm(kk)*ff(kk, rr);
		};

		CHECK(norm(kvec) == Approx(rs.cell().metric().norm(rs.cell().metric().to_covariant(kvec))));
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++) func.hypercubic()[ix][iy][iz][ist] = (ist + 1.0)*gradff(kvec, vec);
				}
			}
		}

		auto div = operations::divergence(func, factor);

		auto diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					for(int ist = 0; ist < func.local_set_size(); ist++){
						auto anvalue = (ist + 1.0)*factor*laplff(kvec, vec);
						diff += fabs(div.hypercubic()[ix][iy][iz][ist] - anvalue);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		CHECK(diff < 1.0e-7);

	}
	
}
#endif

