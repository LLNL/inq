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

    pulay_mixer(const int arg_steps, const double arg_mix_factor):
			size_(-1),
			max_size_(arg_steps),
      mix_factor_(arg_mix_factor){
    }

		template <class mix_type>
    void operator()(mix_type & new_value){

			auto full = (size_ == max_size_);
			size_ = std::min(size_ + 1, max_size_);

			if(size_ == 0){
				//the first step we just store

				ff_ = boost::multi::array<type, 2>({max_size_, new_value.size()});
				dff_= boost::multi::array<type, 2>({max_size_, new_value.size()});
				
				//DATAOPERATIONS LOOP 1D
				for(unsigned ii = 0; ii < new_value.size(); ii++) ff_[0][ii] = new_value[ii];

				return;
			}

			if(size_ == 1){
				//the second step we do linear mixing
				
				//DATAOPERATIONS LOOP 1D
				for(unsigned ii = 0; ii < new_value.size(); ii++){
					dff_[0][ii] = new_value[ii] - ff_[0][ii];
					new_value[ii] = mix_factor_*new_value[ii] + (1.0 - mix_factor_)*ff_[0][ii];
					ff_[0][ii] = new_value[ii];
				}

				return;
			}

			if(full){
				//move all the stored functions, this could be avoided but we do it for simplicity
				for(int istep = 1; istep < size_; istep++){
					//DATAOPERATIONS LOOP 1D
					for(unsigned ii = 0; ii < new_value.size(); ii++){
						ff_[istep - 1][ii] = ff_[istep][ii];
						dff_[istep - 1][ii] = dff_[istep][ii];
					}
				}
			}
			
      //DATAOPERATIONS LOOP 1D
      for(unsigned ii = 0; ii < new_value.size(); ii++){
				ff_[size_ - 1][ii] = new_value[ii];
				dff_[size_ - 1][ii] = new_value[ii] - ff_[size_ - 2][ii];
			}

			boost::multi::array<typename mix_type::value_type, 2> amatrix({size_, size_});

			for(int ii = 0; ii < size_; ii++){
				for(int jj = 0; jj < size_; jj++){
					typename mix_type::value_type aa = 0.0;
					for(unsigned kk = 0; kk < new_value.size(); kk++) aa += conj(dff_[ii][kk])*dff_[jj][kk];
					amatrix[ii][jj] = aa;

					if(ii == jj) std::cout << "norm " << ii << '\t' << aa << std::endl;
					
				}
			}

			//std::cout << amatrix[0][0] << '\t' << amatrix[0][1] << std::endl;
			//std::cout << amatrix[1][0] << '\t' << amatrix[1][1] << std::endl;

			/*
			std::cout << amatrix[0][0] << '\t' << amatrix[0][1] << '\t' << amatrix[0][2] << std::endl;
			std::cout << amatrix[1][0] << '\t' << amatrix[1][1] << '\t' << amatrix[1][2] << std::endl;
			std::cout << amatrix[2][0] << '\t' << amatrix[2][1] << '\t' << amatrix[2][2] << std::endl;
			*/
			
			// REDUCE GRID amatrix

			boost::multi::array<typename mix_type::value_type, 1> alpha(size_, 1.0);

			//std::cout << "alpha = " << alpha[0] << '\t' << alpha[1] << std::endl;
			
			solvers::linear_symmetric(amatrix, alpha);

			//			std::cout << "alpha = " << alpha[0] << '\t' << alpha[1] << std::endl;

			//DATAOPERATIONS
			type sumalpha = 0.0;
			
			//DATAOPERATIONS LOOP 1D
			for(int jj = 0; jj < size_; jj++) sumalpha += alpha[jj];
			//DATAOPERATIONS LOOP 1D
			for(int jj = 0; jj < size_; jj++) alpha[jj] /= sumalpha;

			/*			for(int jj = 0; jj < size_; jj++) alpha[jj] = 0.0;
			alpha[current_] = 0.3;
			alpha[previous_] = 1.0 - 0.3; 
			*/
			
			//DATAOPERATIONS LOOP 1D
      for(unsigned ii = 0; ii < new_value.size(); ii++)	new_value[ii] = 0.0;

			//DATAOPERATIONS LOOP 1D
			for(int jj = 0; jj < size_; jj++) for(unsigned ii = 0; ii < new_value.size(); ii++) new_value[ii] += alpha[jj]*dff_[jj][ii];

			typename mix_type::value_type aa = 0.0;
			for(unsigned kk = 0; kk < new_value.size(); kk++) aa += norm(new_value[kk]);
			std::cout << "norm opt " << aa << std::endl;

			//DATAOPERATIONS LOOP 1D
			for(int jj = 0; jj < size_; jj++) for(unsigned ii = 0; ii < new_value.size(); ii++) new_value[ii] += alpha[jj]*ff_[jj][ii];
			
    }

  private:

		int size_;
		int max_size_;
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
  
  solvers::pulay_mixer<double> lm(5, 0.5);

  v[0] = 0.0;
  v[1] = 22.2;

  lm(v);
  
  REQUIRE(v[0] == 5.0);
  REQUIRE(v[1] == 1.1);
  
}


#endif


#endif
