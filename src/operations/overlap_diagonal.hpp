/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__OVERLAP_DIAGONAL
#define OPERATIONS__OVERLAP_DIAGONAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/array.hpp>
#include <cassert>
#include <operations/integral.hpp>

#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>
#include <gpu/run.hpp>

namespace inq {
namespace operations {

template <class Basis, class PhiMatrix>
gpu::array<typename PhiMatrix::element_type, 1> overlap_diagonal_impl(Basis const & basis, PhiMatrix const & phi1_matrix, PhiMatrix const & phi2_matrix){
	CALI_CXX_MARK_SCOPE("overlap_diagonal(2arg)");

	assert(sizes(phi1_matrix) == sizes(phi2_matrix));

	auto nn = get<1>(phi1_matrix.sizes());

	using type = typename PhiMatrix::element_type;
		
	gpu::array<type, 1> overlap_vector(nn);

	if(nn == 1){
		namespace blas = boost::multi::blas;
		overlap_vector[0] = blas::dot(blas::C(phi1_matrix.rotated()[0]), phi2_matrix.rotated()[0]);
		overlap_vector[0] *= basis.volume_element();
	} else {
		overlap_vector = gpu::run(nn, gpu::reduce(phi1_matrix.size()), zero<type>(),
															[factor = basis.volume_element(), mat1 = begin(phi1_matrix), mat2 = begin(phi2_matrix)] GPU_LAMBDA (auto ist, auto ip) {
																return factor*conj(mat1[ip][ist])*mat2[ip][ist];
															});
	}

	if(basis.comm().size() > 1){
		CALI_CXX_MARK_SCOPE("overlap_diagonal(2arg)::reduce");
		basis.comm().all_reduce_in_place_n(raw_pointer_cast(overlap_vector.data_elements()), overlap_vector.size(), std::plus<>{});
	}
	
	return overlap_vector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type>
auto overlap_diagonal(basis::field_set<Basis, Type> const & phi1, basis::field_set<Basis, Type> const & phi2){
	assert(phi1.basis() == phi2.basis());
	return overlap_diagonal_impl(phi1.basis(), phi1.matrix(), phi2.matrix());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type>
auto overlap_diagonal(states::orbital_set<Basis, Type> const & phi1, states::orbital_set<Basis, Type> const & phi2){
	assert(phi1.basis() == phi2.basis()); 
	return overlap_diagonal_impl(phi1.basis(), phi1.spinor_matrix(), phi2.spinor_matrix());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class field_set_type>
auto overlap_diagonal(const field_set_type & phi){
	CALI_CXX_MARK_SCOPE("overlap_diagonal(1arg)");
	return overlap_diagonal(phi, phi);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
struct value_and_norm {

	value_and_norm() = default;
	
	GPU_FUNCTION value_and_norm(Type const & val, Type const & nrm):
		value(val),
		norm(nrm)
	{
	}

	GPU_FUNCTION value_and_norm(double val):
		value(val),
		norm(val)
	{
	}

	GPU_FUNCTION value_and_norm & operator+=(value_and_norm const & term){
		value += term.value;
		norm  += term.norm;
		return *this;
	}
	
	Type value;
	Type norm;
};

struct identity {
	template <typename Type>
	GPU_FUNCTION auto operator()(Type const tt) const {
		return tt;
	}
};

struct real_part {
	template <typename Type>
	GPU_FUNCTION auto operator()(Type const tt) const {
		return real(tt);
	}
};

template <class Basis, class PhiMatrix, typename Transform>
auto overlap_diagonal_normalized_impl(Basis const & basis, PhiMatrix const & phi1_matrix, PhiMatrix const & phi2_matrix, Transform trans)
	-> gpu::array<decltype(trans(typename PhiMatrix::element_type{})), 1> {

	CALI_CXX_MARK_SCOPE("overlap_diagonal_normalized");

	assert(sizes(phi1_matrix) == sizes(phi2_matrix));

	auto nn = get<1>(phi1_matrix.sizes());

	using type = typename PhiMatrix::element_type;

	auto overlap_and_norm = gpu::run(nn, gpu::reduce(phi1_matrix.size()), zero<value_and_norm<type>>(),
																	 [mat1 = begin(phi1_matrix), mat2 = begin(phi2_matrix)] GPU_LAMBDA (auto ist, auto ip) {
																		 return value_and_norm<type>{conj(mat1[ip][ist])*mat2[ip][ist], conj(mat2[ip][ist])*mat2[ip][ist]};
																	 });
	
	if(basis.comm().size() > 1){
		CALI_CXX_MARK_SCOPE("overlap_diagonal_normalized::reduce");
		basis.comm().all_reduce_in_place_n(reinterpret_cast<type *>(raw_pointer_cast(overlap_and_norm.data_elements())), 2*overlap_and_norm.size(), std::plus<>{});
	}

	gpu::array<decltype(trans(typename PhiMatrix::element_type{})), 1> overlap_vector(nn);

	gpu::run(overlap_vector.size(),
					 [olp = begin(overlap_vector), olpnrm = begin(overlap_and_norm), trans] GPU_LAMBDA (auto ii){
						 olp[ii] = trans(olpnrm[ii].value/olpnrm[ii].norm);
					 });

	return overlap_vector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type, class Transform = identity>
auto overlap_diagonal_normalized(basis::field_set<Basis, Type> const & phi1, basis::field_set<Basis, Type> const & phi2, Transform trans = {}) {
	return overlap_diagonal_normalized_impl(phi1.basis(), phi1.matrix(), phi2.matrix(), trans);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type, class Transform = identity>
auto overlap_diagonal_normalized(states::orbital_set<Basis, Type> const & phi1, states::orbital_set<Basis, Type> const & phi2, Transform trans = {}) {
	return overlap_diagonal_normalized_impl(phi1.basis(), phi1.spinor_matrix(), phi2.spinor_matrix(), trans);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
}
#endif

#ifdef INQ_OPERATIONS_OVERLAP_DIAGONAL_UNIT_TEST
#undef INQ_OPERATIONS_OVERLAP_DIAGONAL_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 800;
	const int nvec = 12;
			
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
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value() + 1);
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1);
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(gpu::array<double, 1>));
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++) {
			auto jjg = aa.set_part().local_to_global(jj).value();
			CHECK(dd[jj] == Approx(-jjg - 1));
		}

		basis::field_set<basis::trivial, double> cc(bas, nvec, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = cc.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value());
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(gpu::array<double, 1>));
								
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) {
				auto jjg = aa.set_part().local_to_global(jj).value();
				CHECK(ee[jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jjg));
			}
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);
			
			CHECK(typeid(decltype(ff)) == typeid(gpu::array<double, 1>));

			CHECK(get<0>(sizes(ff)) == aa.set_part().local_size());
			
			for(int jj = 0; jj < ff.size(); jj++) {
				CHECK(ff[jj] == Approx(dd[jj]/gg[jj]));
			}
 
		}
	}

	SECTION("field_set complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, complex> bb(bas, nvec, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(gpu::array<complex, 1>));

		CHECK(get<0>(sizes(dd)) ==  aa.set_part().local_size());

		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj).value();
			CHECK(fabs(real(dd[jj])) < 1.0e-12);
			CHECK(imag(dd[jj]) == Approx(-jjg - 1));
		}

		basis::field_set<basis::trivial, complex> cc(bas, nvec, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = cc.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(gpu::array<complex, 1>));

			CHECK(get<0>(sizes(ee)) == aa.set_part().local_size());
				
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) {
				auto jjg = aa.set_part().local_to_global(jj).value();
				CHECK(real(ee[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jjg));
			}
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);

			CHECK(typeid(decltype(ff)) == typeid(gpu::array<complex, 1>));

			CHECK(get<0>(sizes(ff)) == aa.set_part().local_size());
				
			for(int jj = 0; jj < ff.size(); jj++) {
				CHECK(fabs(ff[jj]) == Approx(fabs(dd[jj]/gg[jj])));
				CHECK(imag(ff[jj]) == Approx(imag(dd[jj]/gg[jj])));
			}
		}

	}

	SECTION("orbital_set double"){
		
		states::orbital_set<basis::trivial, double> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, double> bb(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value() + 1);
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1);
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(gpu::array<double, 1>));
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++) {
			auto jjg = aa.set_part().local_to_global(jj).value();
			CHECK(dd[jj] == Approx(-jjg - 1));
		}

		states::orbital_set<basis::trivial, double> cc(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = cc.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value());
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(gpu::array<double, 1>));
								
			for(int jj = 0; jj < aa.set_part().local_size(); jj++) {
				auto jjg = aa.set_part().local_to_global(jj).value();
				CHECK(ee[jj] == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jjg));
			}
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);
			
			CHECK(typeid(decltype(ff)) == typeid(gpu::array<double, 1>));

			CHECK(get<0>(sizes(ff)) == aa.set_part().local_size());
			
			for(int jj = 0; jj < ff.size(); jj++) CHECK(ff[jj] == Approx(dd[jj]/gg[jj]));
 
		}
	}

	SECTION("orbital_set complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, complex> bb(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.matrix()[ii][jj] = -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(gpu::array<complex, 1>));
		
		CHECK(get<0>(sizes(dd)) == aa.set_part().local_size());
		
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj).value();
			CHECK(fabs(real(dd[jj])) < 1.0e-12);
			CHECK(imag(dd[jj]) == Approx(-jjg - 1));
		}

		states::orbital_set<basis::trivial, complex> cc(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = cc.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.matrix()[ii][jj] = sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(gpu::array<complex, 1>));

			CHECK(get<0>(sizes(ee)) == aa.set_part().local_size());
				
			for(int jj = 0; jj < ee.size(); jj++) {
				auto jjg = aa.set_part().local_to_global(jj).value();
				CHECK(real(ee[jj]) == Approx(0.5*npoint*(npoint - 1.0)*bas.volume_element()*jjg));
			}
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);

			CHECK(typeid(decltype(ff)) == typeid(gpu::array<complex, 1>));

			CHECK(get<0>(sizes(ff)) == aa.set_part().local_size());

			for(int jj = 0; jj < ff.size(); jj++) {
				CHECK(fabs(ff[jj]) == Approx(fabs(dd[jj]/gg[jj])));
				CHECK(imag(ff[jj]) == Approx(imag(dd[jj]/gg[jj])));
			}
		}

	}
	
	SECTION("orbital_set spinor complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		states::orbital_set<basis::trivial, complex> bb(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);

		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.spinor_set_part().local_size(); jj++){
				auto jjg = aa.spinor_set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.spinor_array()[ii][0][jj] =  20.0*(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				aa.spinor_array()[ii][1][jj] =  -2.0*(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, -M_PI/4 + M_PI/7*iig.value()));
				bb.spinor_array()[ii][0][jj] =  -0.05/(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
				bb.spinor_array()[ii][1][jj] =   0.50/(iig.value() + 1)*sqrt(jjg.value() + 1)*exp(complex(0.0, M_PI/4 + M_PI/7*iig.value()));
			}
		}

		auto dd = operations::overlap_diagonal(aa, bb);
		
		CHECK(typeid(decltype(dd)) == typeid(gpu::array<complex, 1>));

		CHECK(get<0>(sizes(dd)) == aa.spinor_set_part().local_size());
		
		for(int jj = 0; jj < dd.size(); jj++){
			auto jjg = aa.spinor_set_part().local_to_global(jj).value();
			CHECK(fabs(real(dd[jj])) < 1.0e-12);
			CHECK(imag(dd[jj]) == Approx(2.0*(-jjg - 1)));
		}
		
		states::orbital_set<basis::trivial, complex> cc(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.spinor_set_part().local_size(); jj++){
				auto jjg = cc.spinor_set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				cc.spinor_array()[ii][0][jj] = sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));
				cc.spinor_array()[ii][1][jj] = 2.0*sqrt(iig.value())*sqrt(jjg.value())*exp(complex(0.0, M_PI/65.0*iig.value()));        
			}
		}

		{
			auto ee = operations::overlap_diagonal(cc);

			CHECK(typeid(decltype(ee)) == typeid(gpu::array<complex, 1>));

			CHECK(get<0>(sizes(ee)) == aa.spinor_set_part().local_size());
				
			for(int jj = 0; jj < ee.size(); jj++) {
				auto jjg = aa.spinor_set_part().local_to_global(jj).value();
				CHECK(real(ee[jj]) == Approx(2.5*npoint*(npoint - 1.0)*bas.volume_element()*jjg));
			}
		}

		{
			auto ff = operations::overlap_diagonal_normalized(aa, bb);
			auto gg = operations::overlap_diagonal(bb);

			CHECK(typeid(decltype(ff)) == typeid(gpu::array<complex, 1>));

			CHECK(get<0>(sizes(ff)) == aa.spinor_set_part().local_size());

			for(int jj = 0; jj < ff.size(); jj++) {
				CHECK(fabs(ff[jj]) == Approx(fabs(dd[jj]/gg[jj])));
				CHECK(imag(ff[jj]) == Approx(imag(dd[jj]/gg[jj])));
			}
		}
		
	}

}
#endif
