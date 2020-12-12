/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP_DIAGONAL
#define OPERATIONS__OVERLAP_DIAGONAL

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
#include "multi/adaptors/cuda/cublas/context.hpp" // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>
#include <operations/integral.hpp>

#include <caliper/cali.h>

#include <gpu/run.hpp>

namespace inq {
namespace operations {

template <class mat_type>
struct overlap_diagonal_mult {

	mat_type mat1;
	mat_type mat2;
	
	GPU_FUNCTION auto operator()(long ist, long ip) const {
		return conj(mat1[ip][ist])*mat2[ip][ist];
	}
	
};


template <class field_set_type>
math::array<typename field_set_type::element_type, 1> overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

	CALI_CXX_MARK_SCOPE("overlap_diagonal(2arg)");
	
	using type = typename field_set_type::element_type;
		
	math::array<type, 1> overlap_vector(phi1.set_part().local_size());

	assert(size(overlap_vector) == phi2.set_part().local_size());

	if(phi2.set_part().local_size() == 1){
#ifdef ENABLE_CUDA
		if(typeid(typename field_set_type::element_type) == typeid(complex)) {
			cublasZdotc(boost::multi::cuda::cublas::context::get_instance().get(), phi1.basis().part().local_size(),
									(const cuDoubleComplex *) raw_pointer_cast(phi1.matrix().data_elements()), 1, (const cuDoubleComplex *)  raw_pointer_cast(phi2.matrix().data_elements()), 1,
									(cuDoubleComplex *) raw_pointer_cast(overlap_vector.data_elements()));
		} else {
			cublasDdot(boost::multi::cuda::cublas::context::get_instance().get(), phi1.basis().part().local_size(),
								 (const double *) raw_pointer_cast(phi1.matrix().data_elements()), 1, (const double *) raw_pointer_cast(phi2.matrix().data_elements()), 1,
								 (double *) raw_pointer_cast(overlap_vector.data_elements()));
		}
#else

		using boost::multi::blas::dot;
		using boost::multi::blas::conj;
		
		overlap_vector[0] = dot(boost::multi::blas::conj(phi1.matrix().rotated()[0]), phi2.matrix().rotated()[0]);
#endif
		overlap_vector[0] *= phi1.basis().volume_element();
	} else {
		
		overlap_vector = gpu::run(phi1.local_set_size(), gpu::reduce(phi1.basis().part().local_size()), overlap_diagonal_mult<decltype(begin(phi1.matrix()))>{begin(phi1.matrix()), begin(phi2.matrix())});

		gpu::run(phi1.local_set_size(), [ov = begin(overlap_vector), vol_element = phi1.basis().volume_element()] GPU_LAMBDA (auto ist) {
																			ov[ist] = ov[ist]*vol_element;
																		});
	}

	if(phi1.basis().part().parallel()){
		phi1.basis().comm().all_reduce_in_place_n(static_cast<type *>(overlap_vector.data()), overlap_vector.size(), std::plus<>{});
	}
		
	return overlap_vector;
}
	
template <class field_set_type>
auto overlap_diagonal(const field_set_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_diagonal(1arg)");
	
	return overlap_diagonal(phi, phi);
}
	
}
}


#ifdef INQ_OPERATIONS_OVERLAP_DIAGONAL_UNIT_TEST
#undef INQ_OPERATIONS_OVERLAP_DIAGONAL_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::overlap_diagonal", "[operations::overlap_diagonal]") {
	
	using namespace inq;
	using namespace Catch::literals;

	const int npoint = 800;
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
			auto dd = operations::overlap_diagonal(aa, bb);

			CHECK(typeid(decltype(dd)) == typeid(math::array<double, 1>));
				
			for(int jj = 0; jj < nvec; jj++) CHECK(dd[jj] == Approx(-jj));
		}
			
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value());
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa);

			CHECK(typeid(decltype(dd)) == typeid(math::array<double, 1>));
								
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
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value())*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		{
			auto dd = operations::overlap_diagonal(aa, bb);

			CHECK(typeid(decltype(dd)) == typeid(math::array<complex, 1>));

			CHECK(std::get<0>(sizes(dd)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++){
				CHECK(fabs(real(dd[jj])) < 1.0e-12);
				CHECK(imag(dd[jj]) == Approx(-jj));
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
			auto dd = operations::overlap_diagonal(aa);

			CHECK(typeid(decltype(dd)) == typeid(math::array<complex, 1>));
			
			CHECK(std::get<0>(sizes(dd)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++) CHECK(real(dd[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
		}
					
	}

}

#endif
#endif
