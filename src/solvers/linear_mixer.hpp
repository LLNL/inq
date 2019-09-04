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
#include <math/d3vector.hpp>
#include <multi/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/shift.hpp>

namespace solvers {

	template <class mix_type>
  class linear_mixer {

  public:

    linear_mixer(double arg_mix_factor, const mix_type & initial_value):
      mix_factor_(arg_mix_factor),
      old_(initial_value) {
    }

    const auto & operator()(const mix_type & new_value){
      old_ = (1.0 - mix_factor_)*old_ + mix_factor_*new_value;
      return old_;
    }

  private:

    double mix_factor_;
    mix_type old_;
		
	};

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("solvers::linear_mixer", "[solvers::linear_mixer]") {

  solvers::linear_mixer<double> lm(0.5, 10.0);

  REQUIRE(lm(0.0) == 5.0);
  
}


#endif


#endif
