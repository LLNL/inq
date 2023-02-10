/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__ADD
#define OPERATIONS__ADD

/*
 Copyright (C) 2019 Xavier Andrade

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
template <typename FieldType1, typename FieldType2>
void increment(FieldType1 & t1, FieldType2 const & t2){
	assert(t1.basis() == t2.basis());
		
	gpu::run(t1.linear().size(),
					 [t1p = t1.linear().begin(),
						t2p = t2.linear().begin()] GPU_LAMBDA (auto ii){
						 t1p[ii] += t2p[ii];
					 });
}

template <class FVectorType, class BasisType, class Type>
void increment(FVectorType & fvector, basis::field_set<BasisType, Type> const & fset){

	int ist = 0;
	for(auto & field : fvector){
		assert(field.basis() == fset.basis());
		
		gpu::run(fset.basis().local_size(),
						 [fie = begin(field.linear()), fse = begin(fset.matrix()), ist] GPU_LAMBDA (auto ip){
							 fie[ip] += fse[ip][ist];
						 });
		ist++;
	}
}

}
}
#endif

#ifdef INQ_OPERATIONS_ADD_UNIT_TEST
#undef INQ_OPERATIONS_ADD_UNIT_TEST

#include <basis/field.hpp>

#include <catch2/catch_all.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::add", "[operations::add]") {

	using namespace inq;
	using namespace Catch::literals;

	const int N = 100;

	auto comm = boost::mpi3::environment::get_world_instance();
 
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
