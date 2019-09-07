/* -*- indent-tabs-mode: t -*- */

#ifndef SOLVERS_PULAY_MIXER
#define SOLVERS_PULAY_MIXER

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

namespace solvers {

	template <class mix_type>
  class pulay_mixer {

		/*

			The DIIS mixing of Pulay, as described in 
			
			Kresse and Furthmueller, Phys. Rev. B, 54 11169 (1996)

		*/
		
  public:

    pulay_mixer(const int arg_steps, const double arg_mix_factor, const mix_type & initial_value):
			steps_(1),
			max_steps_(arg_steps),
			current_(0),
			previous_(-1),
      mix_factor_(arg_mix_factor),
      ff_({arg_steps, initial_value.size()}),
			dff_({arg_steps, initial_value.size()}) {

      //DATAOPERATIONS
      for(unsigned ii = 0; ii < initial_value.size(); ii++){
				ff_[current_][ii] = initial_value[ii];
				dff_[current_][ii] = initial_value[ii];
			}
    }

    void operator()(mix_type & new_value){

			steps_ = std::min(steps_ + 1, max_steps_);

			current_++;
			if(current_ > steps_) current_ = 0;

			previous_++;
			if(previous_ > steps_) previous_ = 0;
			
      //DATAOPERATIONS
      for(unsigned ii = 0; ii < new_value.size(); ii++){
				ff_[current_][ii] = new_value[ii];
				dff_[current_][ii] = new_value[ii] - ff_[previous_][ii];
			}

    }

  private:

		int steps_;
		int max_steps_;
		int current_;
		int previous_;
    double mix_factor_;
    boost::multi::array<typename mix_type::value_type, 2> ff_;
		boost::multi::array<typename mix_type::value_type, 2> dff_;

		
	};

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("solvers::pulay_mixer", "[solvers::pulay_mixer]") {


  std::vector<double> v(2);

  v[0] =  10.0;
  v[1] = -20.0;
  
  solvers::pulay_mixer<std::vector<double> > lm(5, 0.5, v);

  v[0] = 0.0;
  v[1] = 22.2;

  lm(v);
  
  REQUIRE(v[0] == 5.0);
  REQUIRE(v[1] == 1.1);
  
}


#endif


#endif
