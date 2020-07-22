/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MIXERS__PULAY
#define INQ__MIXERS__PULAY

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
#include <solvers/least_squares.hpp>

namespace inq {
namespace mixers {

template <class Type>
class pulay : public base<Type> {
	
	/*
		The DIIS mixing of Pulay, as described in Kresse and Furthmueller, Phys. Rev. B, 54 11169 (1996)
	*/
		
public:

	pulay(const int arg_steps, const double arg_mix_factor, const long long dim):
		iter_(0),
		max_size_(arg_steps),
		mix_factor_(arg_mix_factor),
		ff_({max_size_, dim}, NAN),
		dff_({max_size_, dim}, NAN){
	}

	void operator()(math::array<Type, 1> & input_value, const math::array<Type, 1> & output_value){

		const double residual_coeff = 0.05;
			
		assert((typename math::array<double, 2>::size_type) input_value.size() == ff_[0].size());
		assert((typename math::array<double, 2>::size_type) output_value.size() == ff_[0].size());

		{
			Type aa = 0.0;
			Type bb = 0.0;
			for(unsigned kk = 0; kk < input_value.size(); kk++){
				aa += fabs(input_value[kk]);
				bb += fabs(output_value[kk]);
			}
			std::cout << "input norm = " << aa << " output norm = " << bb << std::endl;
		}
			
		int size;

		iter_++;

		if(iter_ <= max_size_){
			size = iter_;
		} else {
			size = max_size_;

			//move all the stored functions, this could be avoided but we do it for simplicity
			for(int istep = 1; istep < max_size_; istep++){
				ff_[istep - 1] = ff_[istep];
				dff_[istep - 1] = dff_[istep];
			}
		}

		//DATAOPERATIONS LOOP 1D (one input two outputs)
		for(unsigned ii = 0; ii < input_value.size(); ii++){
			ff_[size - 1][ii] = input_value[ii];
			dff_[size - 1][ii] = output_value[ii] - input_value[ii];
		}
			
		for(int ii = 0; ii < size; ii++){
			Type aa = 0.0;
			Type bb = 0.0;
			for(unsigned kk = 0; kk < input_value.size(); kk++){
				aa += fabs(ff_[ii][kk]);
				bb += conj(dff_[ii][kk])*dff_[ii][kk];
			}
			std::cout << "norm " << ii << '\t' << aa << '\t' << bb << std::endl;
		}
			
		if(iter_ == 1) {
			//DATAOPERATIONS LOOP 1D
			for(unsigned ii = 0; ii < input_value.size(); ii++)	input_value[ii] = (1.0 - mix_factor_)*input_value[ii] + mix_factor_*output_value[ii];

			Type aa = 0.0;
			for(unsigned kk = 0; kk < input_value.size(); kk++) aa += fabs(input_value[kk]);
			std::cout << "norm opt " << aa << std::endl;

			return;
		}

		math::array<Type, 2> amatrix({size + 1, size + 1}, NAN);

		//DATAOPERATIONS LOOP 2D (use overlap)
		for(int ii = 0; ii < size; ii++){
			for(int jj = 0; jj < size; jj++){
				Type aa = 0.0;
				for(unsigned kk = 0; kk < input_value.size(); kk++) aa += conj(dff_[ii][kk])*dff_[jj][kk];
				amatrix[ii][jj] = aa;
			}
		}

		for(int ii = 0; ii < size; ii++){
			amatrix[ii][size] = -1.0;
			amatrix[size][ii] = -1.0;
		}			

		amatrix[size][size] = 0.0;
			
		//std::cout << amatrix[0][0] << '\t' << amatrix[0][1] << std::endl;
		//std::cout << amatrix[1][0] << '\t' << amatrix[1][1] << std::endl;

		/*
			std::cout << amatrix[0][0] << '\t' << amatrix[0][1] << '\t' << amatrix[0][2] << std::endl;
			std::cout << amatrix[1][0] << '\t' << amatrix[1][1] << '\t' << amatrix[1][2] << std::endl;
			std::cout << amatrix[2][0] << '\t' << amatrix[2][1] << '\t' << amatrix[2][2] << std::endl;
		*/
			
		// REDUCE GRID amatrix

		math::array<Type, 1> alpha(size + 1, 0.0);
		alpha[size] = -1.0;

		//std::cout << "alpha = " << alpha[0] << '\t' << alpha[1] << std::endl;
			
		solvers::least_squares(amatrix, alpha);

		//			std::cout << "alpha = " << alpha[0] << '\t' << alpha[1] << std::endl;
			
		//DATAOPERATIONS LOOP
		double sumalpha = 0.0;
		for(int ii = 0; ii < size; ii++) sumalpha += alpha[ii];
			
		//DATAOPERATIONS LOOP
		for(int ii = 0; ii < size; ii++) alpha[ii] /= sumalpha;

		sumalpha = 0.0;
		std::cout << "alpha = ";
		for(int ii = 0; ii < size; ii++) {
			std::cout << alpha[ii] << " "; 
			sumalpha += alpha[ii];
		}
		std::cout << std::endl;
		std::cout << "sumalpha = " << sumalpha << std::endl;

		//			for(int ii = 0; ii < size; ii++) alpha[ii] = 0.0;
		/*			alpha[size - 2] = 1.0 - mix_factor_;
						alpha[size - 1] = mix_factor_;*/
	
		{
			std::fill(input_value.begin(), input_value.end(), 0.0);
				
			for(int jj = 0; jj < size; jj++) for(unsigned ii = 0; ii < input_value.size(); ii++) input_value[ii] += alpha[jj]*dff_[jj][ii];
				
			Type aa = 0.0;
			for(unsigned kk = 0; kk < input_value.size(); kk++) aa += norm(input_value[kk]);
			std::cout << "res norm " << aa << std::endl;
		}

		//DATAOPERATIONS STL FILL
		std::fill(input_value.begin(), input_value.end(), 0.0);
			
		//DATAOPERATIONS LOOP 2D (use gemv)
		for(int jj = 0; jj < size; jj++) {
			for(unsigned ii = 0; ii < input_value.size(); ii++) {
				input_value[ii] += alpha[jj]*(ff_[jj][ii] + residual_coeff*dff_[jj][ii]);
			}
		}

		Type aa = 0.0;
		Type bb = 0.0;
		for(unsigned kk = 0; kk < input_value.size(); kk++) aa += fabs(input_value[kk]);
		for(unsigned kk = 0; kk < input_value.size(); kk++) bb += input_value[kk];
		std::cout << "norm opt " << aa << " " << bb << std::endl;
			
	}

private:

	int iter_;
	int max_size_;
	double mix_factor_;
	math::array<Type, 2> ff_;
	math::array<Type, 2> dff_;
		
};

}
}


#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("mixers::pulay", "[mixers::pulay]") {

	using namespace inq;
	using namespace Catch::literals;
 
  mixers::pulay<double> lm(5, 0.5, 2);

	math::array<double, 1> vin({10.0, -20.0});
	math::array<double, 1> vout({0.0,  22.2});
	
	lm(vin, vout);
  
	CHECK(vin[0] == 5.0_a);
  CHECK(vin[1] == 1.1_a);

	vout = math::array<double, 1>({4.0, 5.5});

	lm(vin, vout);

	CHECK(vin[0] == 4.4216618979_a);
  CHECK(vin[1] == 3.550631855_a);

}


#endif


#endif
