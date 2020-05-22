/* -*- indent-tabs-mode: t -*- */

#ifndef MIXERS__BROYDEN
#define MIXERS__BROYDEN

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
#include <math/vec3d.hpp>
#include <math/array.hpp>
#include <solvers/least_squares.hpp>

namespace mixers {

	template <class Type>
  class broyden : public base<Type> {

  public:

    broyden(const int arg_steps, const double arg_mix_factor, const long long dim):
			iter_(0),
			max_size_(arg_steps),
      mix_factor_(arg_mix_factor),
			dv_({max_size_, dim}, NAN),
			df_({max_size_, dim}, NAN),
			f_old_(dim, NAN),
			vin_old_(dim, NAN),
			last_pos_(-1){
    }

		///////////////////////////////////////////////////////////////////////////////////////////////////
		
		void broyden_extrapolation(math::array<Type, 1> & input_value, int const iter_used, math::array<Type, 1> const & ff){

			double const w0 = 0.01;
			double const ww = 5.0;
						
			if(iter_used == 0){
				for(long ip = 0; ip < input_value.size(); ip++){
					input_value[ip] += mix_factor_*ff[ip];
				}

				return;
			}

			math::array<Type, 2> beta({iter_used, iter_used}, NAN);
			math::array<Type, 1> work(iter_used, NAN);

			for(int ii = 0; ii < iter_used; ii++){
				for(int jj = ii + 1; jj < iter_used; jj++){
					Type aa = 0.0;
					for(unsigned kk = 0; kk < input_value.size(); kk++) aa +=  ww*ww*conj(df_[ii][kk])*df_[jj][kk];
					beta[ii][jj] = aa;
					beta[jj][ii] = conj(aa);
				}
				beta[ii][ii] = w0*w0 + ww*ww;
			}

			for(int ii = 0; ii < iter_used; ii++){
				Type aa = 0.0;
				for(unsigned kk = 0; kk < input_value.size(); kk++) aa += conj(df_[ii][kk])*ff[kk];
				work[ii] = aa;
			}

			//REDUCE beta
			//REDUCE work

			solvers::least_squares(beta, work);
	
			for(long ip = 0; ip < input_value.size(); ip++){
				input_value[ip] += mix_factor_*ff[ip];
			}

			for(int ii = 0; ii < iter_used; ii++){
				for(long ip = 0; ip < input_value.size(); ip++){
					input_value[ip] -= ww*ww*work[ii]*(mix_factor_*df_[ii][ip] + dv_[ii][ip]);
				}
			}
			
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////
		
    void operator()(math::array<Type, 1> & input_value, math::array<Type, 1> const & output_value){

			assert((typename math::array<double, 2>::size_type) input_value.size() == dv_[0].size());
			assert((typename math::array<double, 2>::size_type) output_value.size() == dv_[0].size());

			iter_++;

			math::array<Type, 1> ff(input_value.size());
			
			for(long ip = 0; ip < input_value.size(); ip++){
				ff[ip] = output_value[ip] - input_value[ip]; 
			}

			if(iter_ > 1){

				auto pos = (last_pos_ + 1)%max_size_;
				
				df_[pos] = ff;
				dv_[pos] = input_value;

				for(long ip = 0; ip < input_value.size(); ip++){
					df_[pos][ip] -= f_old_[ip];
					dv_[pos][ip] -= vin_old_[ip];
				}

				gamma_ = 0.0;
				for(long ip = 0; ip < input_value.size(); ip++){
					gamma_ += conj(df_[pos][ip])*df_[pos][ip];
				}
				//REDUCE gamma

				gamma_ = std::max(1e-8, sqrt(gamma_));

				for(long ip = 0; ip < input_value.size(); ip++){
					df_[pos][ip] /= gamma_;
					dv_[pos][ip] /= gamma_;
				}

				last_pos_ = pos;
				
			}

			vin_old_ = input_value;
			f_old_ = ff;

			auto iter_used = std::min(iter_ - 1, max_size_);

			broyden_extrapolation(input_value, iter_used, ff);
				
		}

  private:
		
		int iter_;
		int max_size_;
		double mix_factor_;
		math::array<Type, 2> dv_;
		math::array<Type, 2> df_;
		math::array<Type, 1> f_old_;
		math::array<Type, 1> vin_old_;
		Type gamma_;
		int last_pos_;
		
	};
	
}


#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("mixers::broyden", "[mixers::broyden]") {

	using namespace Catch::literals;
 
  mixers::broyden<double> lm(5, 0.5, 2);

	math::array<double, 1> vin({10.0, -20.0});
	math::array<double, 1> vout({0.0,  22.2});
	
	lm(vin, vout);
  
	CHECK(vin[0] == 5.0_a);
  CHECK(vin[1] == 1.1_a);

	vout = math::array<double, 1>({4.0, 5.5});

	lm(vin, vout);

	CHECK(vin[0] == 4.4419411001_a);
  CHECK(vin[1] == 3.5554591594_a);

}

#endif


#endif
