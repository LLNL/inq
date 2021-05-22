/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP_DIAGONAL
#define OPERATIONS__OVERLAP_DIAGONAL

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
#include <cassert>
#include <operations/integral.hpp>

#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>
#include <gpu/run.hpp>

namespace inq {
namespace operations {

template <class mat_type>
struct overlap_diagonal_mult {

	double factor;
	mat_type mat1;
	mat_type mat2;
	
	GPU_FUNCTION auto operator()(long ist, long ip) const {
		return factor*conj(mat1[ip][ist])*mat2[ip][ist];
	}
	
};


template <class field_set_type>
math::array<typename field_set_type::element_type, 1> overlap_diagonal(const field_set_type & phi1, const field_set_type & phi2){

	CALI_CXX_MARK_SCOPE("overlap_diagonal(2arg)");
	
	using type = typename field_set_type::element_type;
		
	math::array<type, 1> overlap_vector(phi1.set_part().local_size());

	assert(size(overlap_vector) == phi2.set_part().local_size());

	if(phi2.set_part().local_size() == 1){

		namespace blas = boost::multi::blas;
		overlap_vector[0] = blas::dot(blas::C(phi1.matrix().rotated()[0]), phi2.matrix().rotated()[0]);
		overlap_vector[0] *= phi1.basis().volume_element();
	} else {
		
		overlap_vector = gpu::run(phi1.local_set_size(), gpu::reduce(phi1.basis().part().local_size()), overlap_diagonal_mult<decltype(begin(phi1.matrix()))>{phi1.basis().volume_element(), begin(phi1.matrix()), begin(phi2.matrix())});

	}

	if(phi1.basis().comm().size() > 1){
		CALI_CXX_MARK_SCOPE("overlap_diagonal(2arg)::reduce");
		phi1.basis().comm().all_reduce_in_place_n(raw_pointer_cast(overlap_vector.data_elements()), overlap_vector.size(), std::plus<>{});
	}
	
	return overlap_vector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class field_set_type>
auto overlap_diagonal(const field_set_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_diagonal(1arg)");
	
	return overlap_diagonal(phi, phi);
}


template <typename Type>
struct value_and_norm {

	constexpr value_and_norm(Type const & val, Type const & nrm):
		value(val),
		norm(nrm)
	{
	}

	constexpr value_and_norm(double val = 0.0):
		value(val),
		norm(val)
	{
	}

	constexpr auto operator+=(value_and_norm const & term){
		value += term.value;
		norm  += term.norm;
		return *this;
	}
	
	Type value;
	Type norm;
};

template <class mat_type>
struct overlap_diagonal_normalized_mult {

	mat_type mat1;
	mat_type mat2;
	
	GPU_FUNCTION auto operator()(long ist, long ip) const {
		return value_and_norm<decltype(conj(mat1[0][0])*mat2[0][0])>{conj(mat1[ip][ist])*mat2[ip][ist], conj(mat2[ip][ist])*mat2[ip][ist]};
	}
	
};

template <class field_set_type>
math::array<typename field_set_type::element_type, 1> overlap_diagonal_normalized(const field_set_type & phi1, const field_set_type & phi2){

	CALI_CXX_MARK_SCOPE("overlap_diagonal_normalized");

	using type = typename field_set_type::element_type;

	auto overlap_and_norm = gpu::run(phi1.local_set_size(), gpu::reduce(phi1.basis().part().local_size()),
																	 overlap_diagonal_normalized_mult<decltype(begin(phi1.matrix()))>{begin(phi1.matrix()), begin(phi2.matrix())});

	
	if(phi1.basis().comm().size() > 1){
		CALI_CXX_MARK_SCOPE("overlap_diagonal_normalized::reduce");
		phi1.basis().comm().all_reduce_in_place_n(reinterpret_cast<type *>(raw_pointer_cast(overlap_and_norm.data_elements())), 2*overlap_and_norm.size(), std::plus<>{});
	}

	math::array<type, 1> overlap_vector(phi1.set_part().local_size());

	gpu::run(overlap_vector.size(),
					 [olp = begin(overlap_vector), olpnrm = begin(overlap_and_norm)] GPU_LAMBDA (auto ii){
						 olp[ii] = olpnrm[ii].value/olpnrm[ii].norm;
					 });

	return overlap_vector;
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
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value() + 1);
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1);
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(math::array<double, 1>));
		
		for(int jj = 0; jj < nvec; jj++) CHECK(dd[jj] == Approx(-jj - 1));

		basis::field_set<basis::trivial, double> cc(bas, nvec, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = cc.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value());
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(math::array<double, 1>));
								
			for(int jj = 0; jj < nvec; jj++) CHECK(ee[jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);
			
			CHECK(typeid(decltype(ff)) == typeid(math::array<double, 1>));
			
			CHECK(std::get<0>(sizes(ff)) == nvec);
			
			for(int jj = 0; jj < nvec; jj++) CHECK(ff[jj] == Approx(dd[jj]/gg[jj]));
 
		}
	}

	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(math::array<complex, 1>));
		
		CHECK(std::get<0>(sizes(dd)) == nvec);
		
		for(int jj = 0; jj < nvec; jj++){
			CHECK(fabs(real(dd[jj])) < 1.0e-12);
			CHECK(imag(dd[jj]) == Approx(-jj - 1));
		}

		basis::field_set<basis::trivial, complex> cc(bas, nvec, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < nvec; jj++){
				auto jjg = cc.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(math::array<complex, 1>));
			
			CHECK(std::get<0>(sizes(ee)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++) CHECK(real(ee[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jj));
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);

			CHECK(typeid(decltype(ff)) == typeid(math::array<complex, 1>));
			
			CHECK(std::get<0>(sizes(ff)) == nvec);
				
			for(int jj = 0; jj < nvec; jj++) {
				CHECK(fabs(ff[jj]) == Approx(fabs(dd[jj]/gg[jj])));
				CHECK(imag(ff[jj]) == Approx(imag(dd[jj]/gg[jj])));
			}
		}

	}

}

#endif
#endif
