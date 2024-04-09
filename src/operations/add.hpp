/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__ADD
#define OPERATIONS__ADD

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <basis/field.hpp>
#include <basis/field_set.hpp>

#include <cstdlib>
#include <gpu/run.hpp>
#include <algorithm>
#include <functional>
#include <vector>

namespace inq {
namespace operations {

//	Returns a field that has the sum of the values of t1 and t2.
template <class field_type>
field_type add(const field_type & t1, const field_type & t2){
	assert(t1.basis() == t2.basis());
		
	field_type tadd(t1.skeleton());

	gpu::run(t1.linear().size(),
					 [t1p = t1.linear().begin(),
						t2p = t2.linear().begin(),
						taddp = tadd.linear().begin()] GPU_LAMBDA (auto ii){
						 taddp[ii] = t1p[ii] + t2p[ii];
					 });
		
	return tadd;
}
	
//	Returns a field that has the sum of the values of t1, t2 and t3.
template <class field_type>
field_type add(const field_type & t1, const field_type & t2, const field_type & t3){
	assert(t1.basis() == t2.basis());
	assert(t1.basis() == t3.basis());
		
	field_type tadd(t1.skeleton());

	gpu::run(t1.linear().size(),
					 [t1p = t1.linear().begin(),
						t2p = t2.linear().begin(),
						t3p = t3.linear().begin(),
						taddp = tadd.linear().begin()] GPU_LAMBDA (auto ii){
						 taddp[ii] = t1p[ii] + t2p[ii] + t3p[ii];
					 });
	
	return tadd;
}

//	Returns a field_set that has the sum of the values of t1 and t2.
template <typename FieldSetType, typename FieldType>
FieldSetType add(FieldSetType const & t1, FieldType const & t2){
	assert(t1.basis() == t2.basis());
		
	FieldSetType tadd(t1.skeleton());

	gpu::run(t1.local_set_size(), t1.basis().local_size(),
					 [t1p = begin(t1.matrix()), t2p = begin(t2.linear()), taddp = begin(tadd.matrix())] GPU_LAMBDA (auto ist, auto ip){
						 taddp[ip][ist] = t1p[ip][ist] + t2p[ip];
					 });
		
	return tadd;
}

//	Add t2 to t1.
template <typename FieldSet1, typename FieldSet2>
void increment(FieldSet1 & t1, FieldSet2 const & t2){
	assert(t1.basis() == t2.basis());
		
	gpu::run(t1.local_set_size(), t1.basis().local_size(),
					 [t1p = t1.matrix().begin(), t2p = t2.matrix().begin()] GPU_LAMBDA (auto ist, auto ii){
						 t1p[ii][ist] += t2p[ii][ist];
					 });
}

}
}
#endif

#ifdef INQ_OPERATIONS_ADD_UNIT_TEST
#undef INQ_OPERATIONS_ADD_UNIT_TEST

#include <basis/field.hpp>

#include <catch2/catch_all.hpp>
#include <basis/trivial.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

	const int N = 100;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
 
	basis::trivial bas(N, comm);

	SECTION("Add 2 double arrays"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		basis::field<basis::trivial, double> cc(bas);

		aa.fill(1.0);
		bb.fill(2.5);

		cc = operations::add(aa, bb);
		
		for(int ii = 0; ii < cc.linear().size(); ii++) CHECK(cc.linear()[ii] == 3.5_a);

	}
	
	SECTION("Add 2 complex arrays"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		basis::field<basis::trivial, complex> cc(bas);

		aa.fill(complex(1.0, -20.2));
		bb.fill(complex(2.5, 1.2));

		cc = operations::add(aa, bb);
		
		for(int ii = 0; ii < cc.linear().size(); ii++){
			CHECK(real(cc.linear()[ii]) == 3.5_a);
			CHECK(imag(cc.linear()[ii]) == -19.0_a);
		}

	}

	SECTION("Add 3 double arrays"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);
		basis::field<basis::trivial, double> cc(bas);
		basis::field<basis::trivial, double> dd(bas);

		aa.fill(1.0);
		bb.fill(2.5);
		cc.fill(-4.0);

		dd = operations::add(aa, bb, cc);
		
		for(int ii = 0; ii < cc.linear().size(); ii++) CHECK(dd.linear()[ii] == -0.5_a);

	}
	
	SECTION("Add 3 complex arrays"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);
		basis::field<basis::trivial, complex> cc(bas);
		basis::field<basis::trivial, complex> dd(bas);

		aa.fill(complex(1.0, -20.2));
		bb.fill(complex(2.5, 1.2));
		cc.fill(complex(-2.7, 8.6));

		dd = operations::add(aa, bb, cc);
		
		for(int ii = 0; ii < cc.linear().size(); ii++){
			CHECK(real(dd.linear()[ii]) == 0.8_a);
			CHECK(imag(dd.linear()[ii]) == -10.4_a);
		}

	}
	
	SECTION("Increment 2 double arrays"){
		
		basis::field<basis::trivial, double> aa(bas);
		basis::field<basis::trivial, double> bb(bas);

		aa.fill(1.0);
		bb.fill(2.5);

		operations::increment(aa, bb);
		
		for(int ii = 0; ii < aa.linear().size(); ii++) CHECK(aa.linear()[ii] == 3.5_a);

	}
	
	SECTION("Increment 2 complex arrays"){
		
		basis::field<basis::trivial, complex> aa(bas);
		basis::field<basis::trivial, complex> bb(bas);

		aa.fill(complex(1.0, -20.2));
		bb.fill(complex(2.5, 1.2));
		
		operations::increment(aa, bb);
		
		for(int ii = 0; ii < aa.linear().size(); ii++){
			CHECK(real(aa.linear()[ii]) == 3.5_a);
			CHECK(imag(aa.linear()[ii]) == -19.0_a);
		}

	}
	
}
#endif
