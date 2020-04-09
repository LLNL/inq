/* -*- indent-tabs-mode: t -*- */

#ifndef SOLVERS_LINEAR_MIXER
#define SOLVERS_LINEAR_MIXER

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

#include <math/complex.hpp>
#include <math/vec3d.hpp>
#include <math/array.hpp>

namespace solvers {

	template <class type>
  class linear_mixer {

  public:

    linear_mixer(double arg_mix_factor):
      mix_factor_(arg_mix_factor){
    }

		template <class mix_type>
    void operator()(const mix_type & input_value, const mix_type & output_value, mix_type & new_value){
			//note: arguments might alias here			

      //DATAOPERATIONS LOOP 1D
      for(unsigned ii = 0; ii < new_value.size(); ii++){
        new_value[ii] = (1.0 - mix_factor_)*input_value[ii] + mix_factor_*output_value[ii];
      }
    }

  private:
		
    double mix_factor_;
		
	};

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("solvers::linear_mixer", "[solvers::linear_mixer]") {

	using namespace Catch::literals;

  solvers::linear_mixer<double> lm(0.5);

  std::vector<double> vin({10.0, -20.0});
	std::vector<double> vout({0.0,  22.2});
	std::vector<double> vnew(2);  

	lm(vin, vout, vnew);
  
  REQUIRE(vnew[0] == 5.0_a);
  REQUIRE(vnew[1] == 1.1_a);
  
}


#endif


#endif
