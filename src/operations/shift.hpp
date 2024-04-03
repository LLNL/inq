/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SHIFT
#define INQ__OPERATIONS__SHIFT

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/run.hpp>
#include <basis/field_set.hpp>
#include <cassert>

namespace inq {
namespace operations {

template <class array_1d, class FieldSetType1, class FieldSetType2>
void shift(double scale, const array_1d & factor, const FieldSetType1 & shift, FieldSetType2 & phi){

	CALI_CXX_MARK_SCOPE("shift");
	
	assert(size(factor) == phi.spinor_set_part().local_size());

	gpu::run(phi.spinor_set_part().local_size(), phi.spinor_matrix().size(),
					 [factorp = begin(factor), shiftp = begin(shift.spinor_matrix()), phip = begin(phi.spinor_matrix()), scale]
					 GPU_LAMBDA (auto ist, auto ipoint){
						 phip[ipoint][ist] += scale*(factorp[ist]*shiftp[ipoint][ist]);
					 });
}

template <class FieldSetType1, class FieldSetType2>
void shift(typename FieldSetType1::element_type const & factor, const FieldSetType1 & shift, FieldSetType2 & phi){

	CALI_CXX_MARK_SCOPE("shift_const_factor");
	
	//this could be done with axpy
	gpu::run(phi.set_part().local_size(), phi.basis().part().local_size(),
					 [factor, shiftp = begin(shift.matrix()), phip = begin(phi.matrix())]
					 GPU_LAMBDA (auto ist, auto ipoint){
						 phip[ipoint][ist] += factor*shiftp[ipoint][ist];
					 });
}

}
}
#endif

#ifdef INQ_OPERATIONS_SHIFT_UNIT_TEST
#undef INQ_OPERATIONS_SHIFT_UNIT_TEST

#include <basis/trivial.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 185193;
	const int nvec = 7;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	parallel::cartesian_communicator<2> cart_comm(comm, {});

	auto basis_comm = basis::basis_subcomm(cart_comm);
	
	basis::trivial bas(npoint, basis_comm);
	
	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, double> bb(bas, nvec, cart_comm);

		gpu::array<double, 1> factor(aa.set_part().local_size());
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj);
			for(int ii = 0; ii < bas.part().local_size(); ii++){
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 1.0 + 0.765*iig.value()*jjg.value();
				bb.matrix()[ii][jj] = iig.value();
			}
			factor[jj] = 2.0*0.765*jjg.value();
		}

		operations::shift(-0.5, factor, bb, aa);
				
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) CHECK(aa.matrix()[ii][jj] == Approx(1.0));
		}
	}
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		gpu::array<complex, 1> factor(aa.set_part().local_size());
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj);
			for(int ii = 0; ii < bas.part().local_size(); ii++){
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = complex(iig.value(), 1.0 + 0.765*iig.value()*jjg.value());
				bb.matrix()[ii][jj] = iig.value();
			}
			factor[jj] = complex(0.0, 2.0*0.765*jjg.value());
		}

		operations::shift(-0.5, factor, bb, aa);
				
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			auto iig = bas.part().local_to_global(ii);
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) CHECK(real(aa.matrix()[ii][jj]) == Approx(iig.value()));
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) CHECK(imag(aa.matrix()[ii][jj]) == Approx(1.0));
		}
	}	
	
	SECTION("mixed types"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		gpu::array<double, 1> factor(aa.set_part().local_size());
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj);
			for(int ii = 0; ii < bas.part().local_size(); ii++){
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = complex(iig.value(), 1.0 + 0.765*iig.value()*jjg.value());
				bb.matrix()[ii][jj] = complex(0.0, iig.value());
			}
			factor[jj] = 2.0*0.765*jjg.value();
		}

		operations::shift(-0.5, factor, bb, aa);
				
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) {
				auto iig = bas.part().local_to_global(ii);
				CHECK(real(aa.matrix()[ii][jj]) == Approx(iig.value()));
				CHECK(imag(aa.matrix()[ii][jj]) == Approx(1.0));
			}
		}
	}

	SECTION("orbital_set complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, complex> bb(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		
		gpu::array<complex, 1> factor(aa.set_part().local_size());
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj);
			for(int ii = 0; ii < bas.part().local_size(); ii++){
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = complex(iig.value(), 1.0 + 0.765*iig.value()*jjg.value());
				bb.matrix()[ii][jj] = iig.value();
			}
			factor[jj] = complex(0.0, 2.0*0.765*jjg.value());
		}

		operations::shift(-0.5, factor, bb, aa);
				
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			auto iig = bas.part().local_to_global(ii);
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) CHECK(real(aa.matrix()[ii][jj]) == Approx(iig.value()));
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) CHECK(imag(aa.matrix()[ii][jj]) == Approx(1.0));
		}
	}
	
	SECTION("orbital_set spinor complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, complex> bb(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		
		gpu::array<complex, 1> factor(aa.spinor_set_part().local_size());
		
		for(int jj = 0; jj < aa.spinor_set_part().local_size(); jj++){
			auto jjg = aa.spinor_set_part().local_to_global(jj).value();
			for(int ii = 0; ii < bas.part().local_size(); ii++){
				auto iig = bas.part().local_to_global(ii).value();
				aa.spinor_array()[ii][0][jj] = complex(iig, 1.0 + 0.765*iig*jjg);
				bb.spinor_array()[ii][0][jj] = iig;
				aa.spinor_array()[ii][1][jj] = complex(-1.4*iig, 0.2 + 0.765*iig*jjg);
				bb.spinor_array()[ii][1][jj] = iig;
			}
			factor[jj] = complex(0.0, 2.0*0.765*jjg);
		}

		operations::shift(-0.5, factor, bb, aa);
				
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			auto iig = bas.part().local_to_global(ii).value();
			for(int jj = 0; jj < aa.spinor_set_part().local_size(); jj++) {
				CHECK(real(aa.spinor_array()[ii][0][jj]) == Approx(iig));
				CHECK(imag(aa.spinor_array()[ii][0][jj]) == Approx(1.0));
				CHECK(real(aa.spinor_array()[ii][1][jj]) == Approx(-1.4*iig));
				CHECK(imag(aa.spinor_array()[ii][1][jj]) == Approx(0.2));
			}
		}
		
	}	
}
#endif
