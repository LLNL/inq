/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP
#define OPERATIONS__OVERLAP

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/array.hpp>
#include <matrix/gather_scatter.hpp>
#include <operations/integral.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>
#include <parallel/block_array_iterator.hpp>

#include <cassert>

namespace inq {
namespace operations {

template <class Basis, class FullComm, class SetComm, class SetPart1, class MatrixType1, class SetPart2, class MatrixType2>
auto overlap_impl(Basis const & basis, FullComm & full_comm, SetComm & set_comm, SetPart1 const & phi1_set_part, MatrixType1 const & phi1_matrix, SetPart2 const & phi2_set_part, MatrixType2 const & phi2_matrix){

	CALI_CXX_MARK_SCOPE("overlap(2arg)");

	assert(phi1_matrix.size() == phi2_matrix.size());
	
	namespace blas = boost::multi::blas;

	using type = decltype(phi1_matrix[0][0]*phi2_matrix[0][0]);
	
	auto olap = matrix::distributed<type>(full_comm, phi1_set_part.size(), phi2_set_part.size());

	if(olap.comm().size() == 1) {

		olap.block() = blas::gemm(basis.volume_element(), blas::H(phi1_matrix), phi2_matrix);

	} else if(not phi1_set_part.parallel()) {

		auto array = +blas::gemm(basis.volume_element(), blas::H(phi1_matrix), phi2_matrix);
		
		for(int ipart = 0; ipart < olap.partx().comm_size(); ipart++){
			CALI_CXX_MARK_SCOPE("overlap_domains_reduce");

			auto tmp = array({olap.partx().start(ipart), olap.partx().end(ipart)}, {0, olap.party().size()});
			olap.comm().reduce_n(raw_pointer_cast(tmp.base()), tmp.num_elements(), raw_pointer_cast(olap.block().data_elements()), std::plus<>{}, ipart);
		}
		
	} else {

		gpu::array<type, 2> array({phi1_set_part.size(), phi2_set_part.size()}, 0.0);

		for(auto it = parallel::block_array_iterator(phi1_matrix.size(), phi1_set_part, set_comm, phi1_matrix); it != it.end(); ++it){
			auto block = blas::gemm(basis.volume_element(), blas::H(*it), phi2_matrix);
			array({phi1_set_part.start(it.ipart()), phi2_set_part.end(it.ipart())}, {phi1_set_part.start(), phi1_set_part.end()}) = block;
		}
		
		if(full_comm.size() > 1) {
			CALI_CXX_MARK_SCOPE("overlap_reduce");
			full_comm.all_reduce_in_place_n(raw_pointer_cast(array.data_elements()), array.num_elements(), std::plus<>{});
		}

		olap.block() = array({olap.partx().start(), olap.partx().end()}, {olap.party().start(), olap.party().end()});
	}

	return olap;
	
}

template <class Basis, class Type>
auto overlap(basis::field_set<Basis, Type> const & phi1, basis::field_set<Basis, Type> const & phi2){
	assert(phi1.basis() == phi2.basis());
	assert(phi1.full_comm() == phi2.full_comm());
	
	return overlap_impl(phi1.basis(), phi1.full_comm(), phi1.set_comm(), phi1.set_part(), phi1.matrix(), phi2.set_part(), phi2.matrix());
}

template <class Basis, class Type>
auto overlap(states::orbital_set<Basis, Type> const & phi1, states::orbital_set<Basis, Type> const & phi2){
	assert(phi1.basis() == phi2.basis());
	assert(phi1.full_comm() == phi2.full_comm());
	
	return overlap_impl(phi1.basis(), phi1.full_comm(), phi1.set_comm(), phi1.spinor_set_part(), phi1.spinor_matrix(), phi2.spinor_set_part(), phi2.spinor_matrix());
}

template <class FieldSetType>
auto overlap(const FieldSetType & phi){
	CALI_CXX_MARK_SCOPE("overlap(1arg)");
	return overlap(phi, phi);
}



}
}
#endif

#ifdef INQ_OPERATIONS_OVERLAP_UNIT_TEST
#undef INQ_OPERATIONS_OVERLAP_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 100;
	const int nvec = 16;
			
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	auto parstates = comm.size();
	if(comm.size() == 4) parstates = 2;
	if(comm.size() >= 5) parstates = 1;
	
	parallel::cartesian_communicator<2> cart_comm(comm, {boost::mpi3::fill, parstates});
	auto basis_comm = basis::basis_subcomm(cart_comm);

	basis::trivial bas(npoint, basis_comm);

	SECTION("field_set double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, double> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value());
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value());
			}
		}

		{
			auto cc = operations::overlap(aa, bb);
			auto cc_array = matrix::all_gather(cc);
		
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc_array[ii][jj] == Approx(-sqrt(jj)*sqrt(ii)));
			}
		}

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value());
			}
		}

		{
			auto cc = operations::overlap(aa);
			auto cc_array = matrix::all_gather(cc);
			
			CHECK(typeid(decltype(cc_array[0][0])) == typeid(double));

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc_array[ii][jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
			}
		}
		
	}

	SECTION("field_set complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa, bb);
			auto cc_array = matrix::all_gather(cc);         

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) {
					CHECK(fabs(real(cc_array[ii][jj])) < 1.0e-14);
					CHECK(imag(cc_array[ii][jj]) == Approx(sqrt(jj)*sqrt(ii)));
				}
			}
		}

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa);
			auto cc_array = matrix::all_gather(cc);
			
			CHECK(typeid(decltype(cc_array[0][0])) == typeid(complex));

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					CHECK(real(cc_array[ii][jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
					CHECK(fabs(imag(cc_array[ii][jj])) < 1e-13);
				}
			}
		}

	}
	
	SECTION("complex 1x1"){
	
		parallel::cartesian_communicator<2> cart_comm(comm, {comm.size(), 1});
		auto basis_comm = basis::basis_subcomm(cart_comm);
		
		basis::trivial bas(npoint, basis_comm);
	
		const int nvec = 1;
			
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			auto iig = bas.part().local_to_global(ii);
			aa.matrix()[ii][0] = 20.0*sqrt(iig.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
		}
			
		auto cc = operations::overlap(aa);
		auto cc_array = matrix::all_gather(cc);

		using std::get;

		CHECK(get<0>(sizes(cc_array)) == nvec);
		CHECK(get<1>(sizes(cc_array)) == nvec);
			
		CHECK(real(cc_array[0][0]) == Approx(400.0*0.5*npoint*(npoint + 1.0)*bas.volume_element()));
		CHECK(fabs(imag(cc_array[0][0])) < 1e-12);
	}

	SECTION("orbital_set complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, complex> bb(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 2.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.5/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa, bb);
			auto cc_array = matrix::all_gather(cc);

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) {
					CHECK(fabs(real(cc_array[ii][jj])) < 1.0e-14);
					CHECK(imag(cc_array[ii][jj]) == Approx(sqrt(jj)*sqrt(ii)));
				}
			}
		}

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 3.0*sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa);
			auto cc_array = matrix::all_gather(cc);
			
			CHECK(typeid(decltype(cc_array[0][0])) == typeid(complex));

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					CHECK(real(cc_array[ii][jj]) == Approx(4.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
					CHECK(fabs(imag(cc_array[ii][jj])) < 1e-12);
				}
			}
		}

	}

		SECTION("orbital_set spinor complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, complex> bb(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);

		CHECK(aa.spinors());
		CHECK(bb.spinors());        
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
				auto jjg = aa.spinor_set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.spinor_array()[ii][0][jj] = 2.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
				aa.spinor_array()[ii][1][jj] = 10.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
				bb.spinor_array()[ii][0][jj] = -0.5/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.spinor_array()[ii][1][jj] = -0.1/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));                
			}
		}
		
		{
			auto cc = operations::overlap(aa, bb);
			auto cc_array = matrix::all_gather(cc);

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) {
					CHECK(fabs(real(cc_array[ii][jj])) < 1.0e-12);
					CHECK(imag(cc_array[ii][jj]) == Approx(2.0*sqrt(jj)*sqrt(ii)));
				}
			}
		}

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
				auto jjg = aa.spinor_set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.spinor_array()[ii][0][jj] = 3.0*sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
				aa.spinor_array()[ii][1][jj] = 2.0*sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa);
			auto cc_array = matrix::all_gather(cc);
			
			CHECK(typeid(decltype(cc_array[0][0])) == typeid(complex));

			using std::get;

			CHECK(get<0>(sizes(cc_array)) == nvec);
			CHECK(get<1>(sizes(cc_array)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					CHECK(real(cc_array[ii][jj]) == Approx(6.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
					CHECK(fabs(imag(cc_array[ii][jj])) < 2e-12);
				}
			}
		}
		
	}

}
#endif
