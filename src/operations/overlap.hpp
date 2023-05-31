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

#include <cassert>

namespace inq {
namespace operations {

template <class FieldSetType1, class FieldSetType2>
auto overlap(const FieldSetType1 & phi1, const FieldSetType2 & phi2){

	CALI_CXX_MARK_SCOPE("overlap(2arg)");

	namespace blas = boost::multi::blas;

	auto matrix = matrix::distributed<typename FieldSetType1::element_type>(phi1.full_comm(), phi2.set_size(), phi1.set_size());

	if(matrix.comm().size() == 1) {

		matrix.block() = blas::gemm(phi1.basis().volume_element(), blas::H(phi2.matrix()), phi1.matrix());

	} else if(not phi1.set_part().parallel()) {

		auto array = +blas::gemm(phi1.basis().volume_element(), blas::H(phi2.matrix()), phi1.matrix());
		
		for(int ipart = 0; ipart < matrix.partx().comm_size(); ipart++){
			CALI_CXX_MARK_SCOPE("overlap_domains_reduce");

			auto tmp = array({matrix.partx().start(ipart), matrix.partx().end(ipart)}, {0, matrix.party().size()});
			matrix.comm().reduce_n(raw_pointer_cast(tmp.base()), tmp.num_elements(), raw_pointer_cast(matrix.block().data_elements()), std::plus<>{}, ipart);
		}
		
	} else {

		gpu::array<typename FieldSetType1::element_type, 2> array({phi2.set_size(), phi1.set_size()}, 0.0);

		for(auto it = phi1.par_set_begin(); it != phi1.par_set_end(); ++it){
			auto block = blas::gemm(phi1.basis().volume_element(), blas::H(phi2.matrix()), it.matrix());
			array({phi2.set_part().start(), phi2.set_part().end()}, {phi1.set_part().start(it.set_ipart()), phi1.set_part().end(it.set_ipart())}) = block;
		}
		
		if(phi1.full_comm().size() > 1) {
			CALI_CXX_MARK_SCOPE("overlap_reduce");
			phi1.full_comm().all_reduce_in_place_n(raw_pointer_cast(array.data_elements()), array.num_elements(), std::plus<>{});
		}

		matrix::scatter(array, matrix, /* root = */ 0);		
	}

	return matrix;
	
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
	if(comm.size() >= 5) parstates = 1;
	
	parallel::cartesian_communicator<2> cart_comm(comm, {boost::mpi3::fill, parstates});
	auto basis_comm = basis::basis_subcomm(cart_comm);

	basis::trivial bas(npoint, basis_comm);

	SECTION("double"){
		
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
			CHECK(std::get<0>(sizes(cc_array)) == nvec);
			CHECK(std::get<1>(sizes(cc_array)) == nvec);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc_array[ii][jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
			}
		}
		
	}

	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.local_set_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa, bb);
			auto cc_array = matrix::all_gather(cc);			

			CHECK(std::get<0>(sizes(cc_array)) == nvec);
			CHECK(std::get<1>(sizes(cc_array)) == nvec);
				
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
			CHECK(std::get<0>(sizes(cc_array)) == nvec);
			CHECK(std::get<1>(sizes(cc_array)) == nvec);
				
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
			
		CHECK(std::get<0>(sizes(cc_array)) == nvec);
		CHECK(std::get<1>(sizes(cc_array)) == nvec);
			
		CHECK(real(cc_array[0][0]) == Approx(400.0*0.5*npoint*(npoint + 1.0)*bas.volume_element()));
		CHECK(fabs(imag(cc_array[0][0])) < 1e-12);
	}

}
#endif
