/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MIXERS__LINEAR
#define INQ__MIXERS__LINEAR

/*
 Copyright (C) 2019-2020 Xavier Andrade

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

#include <math/complex.hpp>
#include <math/vector3.hpp>
#include <math/array.hpp>
#include <mixers/base.hpp>

namespace inq {
namespace mixers {

template <class Type>
class linear : public base<Type> {

public:
	
	linear(double arg_mix_factor):
		mix_factor_(arg_mix_factor){
	}
	
	void operator()(math::array<Type, 1> & input_value, math::array<Type, 1>  const & output_value){
		//note: arguments might alias here
		
		CALI_CXX_MARK_SCOPE("linear_mixing");
		
		gpu::run(input_value.size(),
						 [iv = begin(input_value), ov = begin(output_value), mix = mix_factor_] GPU_LAMBDA (auto ii){
							 iv[ii] = (1.0 - mix)*iv[ii] + mix*ov[ii];
						 });
	}

private:
		
	double mix_factor_;
		
};

}
}

#ifdef INQ_MIXERS_LINEAR_UNIT_TEST
#undef INQ_MIXERS_LINEAR_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("mixers::linear", "[mixers::linear]") {

	using namespace inq;
	using namespace Catch::literals;

  mixers::linear<double> lm(0.5);

  math::array<double, 1> vin({10.0, -20.0});
	math::array<double, 1> vout({0.0,  22.2});

	lm(vin, vout);
  
  CHECK(vin[0] == 5.0_a);
  CHECK(vin[1] == 1.1_a);
  
}


#endif


#endif
