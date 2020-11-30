/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP
#define OPERATIONS__OVERLAP

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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
#include <cassert>
#ifdef ENABLE_CUDA
#include "multi/adaptors/blas/cuda.hpp" // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>
#include <operations/integral.hpp>

#include <caliper/cali.h>

namespace inq {
namespace operations {

template <class field_set_type>
auto overlap(const field_set_type & phi1, const field_set_type & phi2){

	CALI_CXX_MARK_SCOPE("overlap(2arg)");
	
	// no state parallelization for now
	assert(not phi1.set_part().parallel());
		
	using boost::multi::blas::gemm;
	using boost::multi::blas::hermitized;

	auto overlap_matrix = gemm(phi1.basis().volume_element(), hermitized(phi2.matrix()), phi1.matrix());

	if(phi1.basis().part().parallel()){
		phi1.basis().comm().all_reduce_in_place_n(static_cast<typename field_set_type::element_type *>(overlap_matrix.data()), overlap_matrix.num_elements(), std::plus<>{});
	}
		
	return overlap_matrix;
}

template <class field_set_type>
auto overlap(const field_set_type & phi){

	CALI_CXX_MARK_SCOPE("overlap(1arg)");
 
	// no state parallelization for now
	assert(not phi.set_part().parallel());

	using boost::multi::blas::herk;
	using boost::multi::blas::hermitized;
		
	auto overlap_matrix = herk(phi.basis().volume_element(), hermitized(phi.matrix()));

	if(phi.basis().part().parallel()){
		phi.basis().comm().all_reduce_in_place_n(static_cast<typename field_set_type::element_type *>(overlap_matrix.data()), overlap_matrix.num_elements(), std::plus<>{});
	}
		
	return overlap_matrix;
}

template <class field_type>
auto overlap_single(const field_type & phi1, const field_type & phi2){
	CALI_CXX_MARK_SCOPE("overlap_single(2arg)");
	
	return integral(phi1, phi2, [](auto t1, auto t2){ return conj(t1)*t2; });
}
	
template <class field_type>
auto overlap_single(field_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_single(1arg)");
			
	return overlap_single(phi, phi);
}
	
}
}


#ifdef INQ_OPERATIONS_OVERLAP_UNIT_TEST
#undef INQ_OPERATIONS_OVERLAP_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::overlap", "[operations::overlap]") {
	
	using namespace inq;
	using namespace Catch::literals;

	const int npoint = 100;
	const int nvec = 12;
			
	auto comm = boost::mpi3::environment::get_world_instance();
		
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {1, comm.size()});

	auto basis_comm = cart_comm.axis(1);

	CHECK(basis_comm.size() == comm.size());
		
	basis::trivial bas(npoint, basis_comm);

	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, double> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value());
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value());
			}
		}

		{
			auto cc = operations::overlap(aa, bb);

			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc[ii][jj] == Approx(-sqrt(jj)*sqrt(ii)));
			}
		}

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value());
			}
		}

		/* This is disabled because it causes problems with nvcc and multi.
				 
			 {
			 auto cc = operations::overlap(aa);

			 CHECK(std::get<0>(sizes(cc)) == nvec);
			 CHECK(std::get<1>(sizes(cc)) == nvec);
				
			 for(int ii = 0; ii < nvec; ii++){
			 for(int jj = 0; jj < nvec; jj++) CHECK(cc[ii][jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
			 }
			 }
		*/

	}

	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa, bb);

			CHECK(std::get<0>(sizes(cc)) == nvec);
			CHECK(std::get<1>(sizes(cc)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) {
					CHECK(fabs(real(cc[ii][jj])) < 1.0e-14);
					CHECK(imag(cc[ii][jj]) == Approx(sqrt(jj)*sqrt(ii)));
				}
			}
		}

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto cc = operations::overlap(aa);

			CHECK(std::get<0>(sizes(cc)) == nvec);
			CHECK(std::get<1>(sizes(cc)) == nvec);
				
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					CHECK(real(cc[ii][jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*sqrt(jj)*sqrt(ii)) );
					CHECK(fabs(imag(cc[ii][jj])) < 1e-13);
				}
			}
		}

	}

		
	SECTION("Overlap single double"){
			
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
			
		aa = 2.0;
		bb = 0.8;
		
		CHECK(operations::overlap_single(aa, bb) == 1.6_a);
			
		for(int ii = 0; ii < bas.part().local_size(); ii++)	{
			auto iig = bas.part().local_to_global(ii);
			aa.linear()[ii] = pow(iig.value() + 1, 2);
			bb.linear()[ii] = 1.0/(iig.value() + 1);
		}
			
		CHECK(operations::overlap_single(aa, bb) == Approx(0.5*npoint*(npoint + 1.0)*bas.volume_element()));
			
	}
		
	SECTION("Integral product complex"){
			
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
			
		aa = complex(2.0, -0.3);
		bb = complex(0.8, 0.01);
		
		CHECK(real(operations::overlap_single(aa, bb)) == 1.597_a);
		CHECK(imag(operations::overlap_single(aa, bb)) == 0.26_a);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++)	{
			auto iig = bas.part().local_to_global(ii);
			aa.linear()[ii] = pow(iig.value() + 1, 2)*exp(complex(0.0, -M_PI/8 + 2.0*M_PI/(iig.value() + 1)));
			bb.linear()[ii] = 1.0/(iig.value() + 1)*exp(complex(0.0, M_PI/8 + 2.0*M_PI/(iig.value() + 1)));
		}

		CHECK(real(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*npoint*(npoint + 1.0)*bas.volume_element()));
		CHECK(imag(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*npoint*(npoint + 1.0)*bas.volume_element()));
		
	}

	SECTION("complex 1x1"){

		const int nvec = 1;
			
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			auto iig = bas.part().local_to_global(ii);
			aa.matrix()[ii][0] = 20.0*sqrt(iig.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
		}
			
		auto cc = operations::overlap(aa);
			
		CHECK(std::get<0>(sizes(cc)) == nvec);
		CHECK(std::get<1>(sizes(cc)) == nvec);
			
		CHECK(real(cc[0][0]) == Approx(400.0*0.5*npoint*(npoint + 1.0)*bas.volume_element()));
		CHECK(fabs(imag(cc[0][0])) < 1e-13);

	}

}

#endif
#endif
