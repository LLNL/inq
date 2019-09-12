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
#include <solvers/linear.hpp>

namespace solvers {

	template <class type>
  class pulay_mixer {

		/*

			The DIIS mixing of Pulay, as described in Kresse and Furthmueller, Phys. Rev. B, 54 11169 (1996)

		*/
		
  public:

		template <class mix_type>
    pulay_mixer(const int arg_steps, const double arg_mix_factor, const mix_type & initial_value):
			size_(0),
			max_size_(arg_steps),
			current_(-1),
			previous_(-2),
      mix_factor_(arg_mix_factor),
      ff_({arg_steps, initial_value.size()}),
			dff_({arg_steps, initial_value.size()}) {

      //DATAOPERATIONS
      for(unsigned ii = 0; ii < initial_value.size(); ii++)	ff_[0][ii] = initial_value[ii];
    }

		template <class mix_type>
    void operator()(mix_type & new_value){
			
			size_ = std::min(size_ + 1, max_size_);

			current_++;
			if(current_ >= max_size_) current_ = 0;

			previous_++;
			if(previous_ >= max_size_) previous_ = 0;

			if(size_ == 1){
				//the first step we do linear mixing
				
				//DATAOPERATIONS
				for(unsigned ii = 0; ii < new_value.size(); ii++){
					dff_[0][ii] = new_value[ii] - ff_[0][ii];
					new_value[ii] = mix_factor_*new_value[ii] + (1.0 - mix_factor_)*ff_[0][ii];
					ff_[0][ii] = new_value[ii];
				}

				return;
			}

			//std::cout << "size " << size_ << '\t' << current_ << '\t' << previous_<< std::endl;

      //DATAOPERATIONS
      for(unsigned ii = 0; ii < new_value.size(); ii++){
				ff_[current_][ii] = new_value[ii];
				dff_[current_][ii] = new_value[ii] - ff_[previous_][ii];
			}

			boost::multi::array<typename mix_type::value_type, 2> amatrix({size_, size_});

			for(int ii = 0; ii < size_; ii++){
				for(int jj = 0; jj < size_; jj++){
					typename mix_type::value_type aa = 0.0;
					for(unsigned kk = 0; kk < new_value.size(); kk++) aa += conj(dff_[ii][kk])*dff_[jj][kk];
					amatrix[ii][jj] = aa;
				}
			}

			//std::cout << amatrix[0][0] << '\t' << amatrix[0][1] << std::endl;
			//std::cout << amatrix[1][0] << '\t' << amatrix[1][1] << std::endl;
			
			// REDUCE GRID amatrix

			boost::multi::array<typename mix_type::value_type, 1> alpha(size_, 1.0);

			//std::cout << "alpha = " << alpha[0] << '\t' << alpha[1] << std::endl;
			
			solvers::linear_symmetric(amatrix, alpha);

			//			std::cout << "alpha = " << alpha[0] << '\t' << alpha[1] << std::endl;

			//DATAOPERATIONS
			type sumalpha = 0.0;
			for(int jj = 0; jj < size_; jj++) sumalpha += alpha[jj];
			for(int jj = 0; jj < size_; jj++) alpha[jj] /= sumalpha;
		
			//DATAOPERATIONS
      for(unsigned ii = 0; ii < new_value.size(); ii++)	new_value[ii] = 0.0;
			for(int jj = 0; jj < size_; jj++) for(unsigned ii = 0; ii < new_value.size(); ii++) new_value[ii] += alpha[jj]*ff_[jj][ii];
			
    }

  private:

		int size_;
		int max_size_;
		int current_;
		int previous_;
    double mix_factor_;
    boost::multi::array<type, 2> ff_;
		boost::multi::array<type, 2> dff_;
		
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
  
  solvers::pulay_mixer<double> lm(5, 0.5, v);

  v[0] = 0.0;
  v[1] = 22.2;

  lm(v);
  
  REQUIRE(v[0] == 5.0);
  REQUIRE(v[1] == 1.1);
  
}


#endif


#endif
