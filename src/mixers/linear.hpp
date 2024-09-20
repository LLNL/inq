/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MIXERS__LINEAR
#define INQ__MIXERS__LINEAR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/complex.hpp>
#include <math/vector3.hpp>
#include <gpu/array.hpp>
#include <mixers/base.hpp>

namespace inq {
namespace mixers {

template <class FieldType>
class linear : public base<FieldType> {

public:
	
	linear(double arg_mix_factor):
		mix_factor_(arg_mix_factor){
	}
	
	void operator()(FieldType & input_value, FieldType const & output_value){
		//note: arguments might alias here

		assert(input_value.basis() == output_value.basis());
		assert(input_value.local_set_size() == output_value.local_set_size());
		
		CALI_CXX_MARK_SCOPE("linear_mixing");
		
		gpu::run(input_value.local_set_size(), input_value.basis().size(),
						 [iv = begin(input_value.matrix()), ov = begin(output_value.matrix()), mix = mix_factor_] GPU_LAMBDA (auto is, auto ip){
							 iv[ip][is] = (1.0 - mix)*iv[ip][is] + mix*ov[ip][is];
						 });

	}

private:
		
	double mix_factor_;
		
};

}
}
#endif

#ifdef INQ_MIXERS_LINEAR_UNIT_TEST
#undef INQ_MIXERS_LINEAR_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

	basis::trivial bas(2, parallel::communicator{boost::mpi3::environment::get_self_instance()});
	
	basis::field_set<basis::trivial, double> vin(bas, 2);
	vin.matrix()[0][0] = 10.0;
	vin.matrix()[1][0] = -20.0;
	vin.matrix()[0][1] = 3.45;
	vin.matrix()[1][1] = 192.34;

	basis::field_set<basis::trivial, double> vout(bas, 2);
	vout.matrix()[0][0] = 0.0;
	vout.matrix()[1][0] = 22.2;
	vout.matrix()[0][1] = 1.87;
	vout.matrix()[1][1] = -133.55;

	mixers::linear<decltype(vin)> mixer(0.5);
	
	mixer(vin, vout);
  
	CHECK(vin.matrix()[0][0] == 5.0_a);
	CHECK(vin.matrix()[1][0] == 1.1_a);
	CHECK(vin.matrix()[0][1] == 2.66_a);
	CHECK(vin.matrix()[1][1] == 29.395_a);
  
}
#endif
