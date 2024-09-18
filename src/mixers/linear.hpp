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

template <class ArrayType>
class linear : public base<ArrayType> {

public:
	
	linear(double arg_mix_factor):
		mix_factor_(arg_mix_factor){
	}
	
	void operator()(ArrayType & input_value, ArrayType const & output_value){
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

  gpu::array<double, 1> vin({10.0, -20.0});
	gpu::array<double, 1> vout({0.0,  22.2});

  mixers::linear<decltype(vin)> lm(0.5);
	
	lm(vin, vout);
  
  CHECK(vin[0] == 5.0_a);
  CHECK(vin[1] == 1.1_a);
  
}
#endif
