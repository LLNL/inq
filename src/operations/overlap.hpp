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

	CALI_CXX_MARK_SCOPE("overlap 2 arg");
	
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

	CALI_CXX_MARK_SCOPE("overlap 1 arg");
 
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

template <class field_set_type>
math::array<typename field_set_type::element_type, 1> overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

	CALI_CXX_MARK_SCOPE("overlap_diagonal 2 arg");
	
	using type = typename field_set_type::element_type;
		
	math::array<type, 1> overlap_vector(phi1.set_part().local_size());

	assert(size(overlap_vector) == phi2.set_part().local_size());

	if(phi2.set_part().local_size() == 1){
#ifdef ENABLE_CUDA
		cublasHandle_t handle;
		cublasCreate(&handle);
		if(typeid(typename field_set_type::element_type) == typeid(complex)) {
			cublasZdotc(handle, phi1.basis().part().local_size(),
									(const cuDoubleComplex *) raw_pointer_cast(phi1.matrix().data_elements()), 1, (const cuDoubleComplex *)  raw_pointer_cast(phi2.matrix().data_elements()), 1,
									(cuDoubleComplex *) raw_pointer_cast(overlap_vector.data_elements()));
		} else {
			cublasDdot(handle, phi1.basis().part().local_size(),
								 (const double *) raw_pointer_cast(phi1.matrix().data_elements()), 1, (const double *) raw_pointer_cast(phi2.matrix().data_elements()), 1,
								 (double *) raw_pointer_cast(overlap_vector.data_elements()));
		}
		cublasDestroy(handle);
#else

		using boost::multi::blas::dot;
		using boost::multi::blas::conj;
		
		overlap_vector[0] = dot(boost::multi::blas::conj(phi1.matrix().rotated()[0]), phi2.matrix().rotated()[0]);
#endif
		overlap_vector[0] *= phi1.basis().volume_element();
	} else {
		
		//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifndef ENABLE_CUDA
		
		//OPTIMIZATION: this can be done more efficiently
		for(int ii = 0; ii < phi1.local_set_size(); ii++){
			type aa = 0.0;
			for(int ip = 0; ip < phi1.basis().part().local_size(); ip++) aa += conj(phi1.matrix()[ip][ii])*phi2.matrix()[ip][ii];
			overlap_vector[ii] = aa*phi1.basis().volume_element();
		}
		
#else
		
		{
			auto npoints = phi1.basis().part().local_size();
			auto vol_element = phi1.basis().volume_element();
			auto phi1p = begin(phi1.matrix());
			auto phi2p = begin(phi2.matrix());
			auto overlap = begin(overlap_vector);
			
			//OPTIMIZATION: here we should parallelize over points as well 
			gpu::run(phi1.local_set_size(),
							 [=] __device__ (auto ist){
								 type aa = 0.0;
								 for(int ip = 0; ip < npoints; ip++){
									 auto p1 = phi1p[ip][ist];
									 auto p2 = phi2p[ip][ist];
									 aa += conj(p1)*p2;
									 
								 }
								 
								 overlap[ist] = vol_element*aa;
							 });
		}
		
#endif

	}
	
	if(phi1.basis().part().parallel()){
		phi1.basis().comm().all_reduce_in_place_n(static_cast<type *>(overlap_vector.data()), overlap_vector.size(), std::plus<>{});
	}
		
	return overlap_vector;
}
	
template <class field_set_type>
auto overlap_diagonal(const field_set_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_diagonal 1 arg");
	
	return overlap_diagonal(phi, phi);
}
	
template <class field_type>
auto overlap_single(const field_type & phi1, const field_type & phi2){
	CALI_CXX_MARK_SCOPE("overlap_single 2 arg");
	
	return integral(phi1, phi2, [](auto t1, auto t2){ return conj(t1)*t2; });
}
	
template <class field_type>
auto overlap_single(field_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_single 1 arg");
			
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
				aa.matrix()[ii][jj] = 20.0*(iig + 1)*sqrt(jjg);
				bb.matrix()[ii][jj] = -0.05/(iig + 1)*sqrt(jjg);
			}
		}

		{
			auto cc = operations::overlap(aa, bb);

			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++) CHECK(cc[ii][jj] == Approx(-sqrt(jj)*sqrt(ii)));
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa, bb);
				
			for(int jj = 0; jj < nvec; jj++) CHECK(dd[jj] == Approx(-jj));
		}
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig)*sqrt(jjg);
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

		{
			auto dd = operations::overlap_diagonal(aa);
								
			for(int jj = 0; jj < nvec; jj++) CHECK(dd[jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
		}
					
			
	}

	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig + 1)*sqrt(jjg)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig));
				bb.matrix()[ii][jj] = -0.05/(iig + 1)*sqrt(jjg)*exp(complex(0.0, M_PI/4 + M_PI/7*iig));
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

		{
			auto dd = operations::overlap_diagonal(aa, bb);

			CHECK(std::get<0>(sizes(dd)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++){
				CHECK(fabs(real(dd[jj])) < 1.0e-14);
				CHECK(imag(dd[jj]) == Approx(-jj));
			}
		}
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig)*sqrt(jjg)*exp(complex(0.0, M_PI/65.0*iig));
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

		{
			auto dd = operations::overlap_diagonal(aa);

			CHECK(std::get<0>(sizes(dd)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++) CHECK(real(dd[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
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
			aa.linear()[ii] = pow(iig + 1, 2);
			bb.linear()[ii] = 1.0/(iig + 1);
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
			aa.linear()[ii] = pow(iig + 1, 2)*exp(complex(0.0, -M_PI/8 + 2.0*M_PI/(iig + 1)));
			bb.linear()[ii] = 1.0/(iig + 1)*exp(complex(0.0, M_PI/8 + 2.0*M_PI/(iig + 1)));
		}

		CHECK(real(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*npoint*(npoint + 1.0)*bas.volume_element()));
		CHECK(imag(operations::overlap_single(aa, bb)) == Approx(sqrt(2.0)*0.25*npoint*(npoint + 1.0)*bas.volume_element()));
		
	}

	SECTION("complex 1x1"){

		const int nvec = 1;
			
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			auto iig = bas.part().local_to_global(ii);
			aa.matrix()[ii][0] = 20.0*sqrt(iig + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig));
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
