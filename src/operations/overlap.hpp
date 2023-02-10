/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP
#define OPERATIONS__OVERLAP

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa.

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

#include <inq_config.h>

#include <math/array.hpp>
#include <math/subspace_matrix.hpp>
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

	if(not phi1.set_part().parallel()){
		
		auto overlap_matrix =+ blas::gemm(phi1.basis().volume_element(), blas::H(phi2.matrix()), phi1.matrix());
		
		if(phi1.basis().comm().size() > 1) {
			CALI_CXX_MARK_SCOPE("overlap(2arg)_mpi_reduce");	
			phi1.basis().comm().all_reduce_in_place_n(raw_pointer_cast(overlap_matrix.data_elements()), overlap_matrix.num_elements(), std::plus<>{});
		}

		return math::subspace_matrix<typename decltype(overlap_matrix)::element_type>(phi1.full_comm(), std::move(overlap_matrix));
		
	} else {

		math::subspace_matrix<typename FieldSetType1::element_type> overlap_matrix(phi1.full_comm(), phi1.set_size(), 0.0);

		for(auto it = phi1.par_set_begin(); it != phi1.par_set_end(); ++it){
			auto block = blas::gemm(phi1.basis().volume_element(), blas::H(phi2.matrix()), it.matrix());
			overlap_matrix.array()({phi2.set_part().start(), phi2.set_part().end()}, {phi1.set_part().start(it.set_ipart()), phi1.set_part().end(it.set_ipart())}) = block;
		}
		
		if(phi1.full_comm().size() > 1) {
			CALI_CXX_MARK_SCOPE("overlap(2arg)_mpi_reduce");	
			phi1.full_comm().all_reduce_in_place_n(raw_pointer_cast(overlap_matrix.array().data_elements()), overlap_matrix.array().num_elements(), std::plus<>{});
		}
		
		return overlap_matrix;
	}
	

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
			
	auto comm = boost::mpi3::environment::get_world_instance();

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

			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc.array()[ii][jj] == Approx(-sqrt(jj)*sqrt(ii)));
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

			CHECK(typeid(decltype(cc.array()[0][0])) == typeid(double));
			CHECK(std::get<0>(sizes(cc.array())) == nvec);
			CHECK(std::get<1>(sizes(cc.array())) == nvec);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc.array()[ii][jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
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
			

			CHECK(std::get<0>(sizes(cc.array())) == nvec);
			CHECK(std::get<1>(sizes(cc.array())) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) {
					CHECK(fabs(real(cc.array()[ii][jj])) < 1.0e-14);
					CHECK(imag(cc.array()[ii][jj]) == Approx(sqrt(jj)*sqrt(ii)));
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

			CHECK(typeid(decltype(cc.array()[0][0])) == typeid(complex));
			CHECK(std::get<0>(sizes(cc.array())) == nvec);
			CHECK(std::get<1>(sizes(cc.array())) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					CHECK(real(cc.array()[ii][jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
					CHECK(fabs(imag(cc.array()[ii][jj])) < 1e-13);
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
			
		CHECK(std::get<0>(sizes(cc.array())) == nvec);
		CHECK(std::get<1>(sizes(cc.array())) == nvec);
			
		CHECK(real(cc.array()[0][0]) == Approx(400.0*0.5*npoint*(npoint + 1.0)*bas.volume_element()));
		CHECK(fabs(imag(cc.array()[0][0])) < 1e-13);

		}

}
#endif
