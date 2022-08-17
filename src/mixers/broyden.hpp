/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MIXERS__BROYDEN
#define INQ__MIXERS__BROYDEN

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
#include <mixers/base.hpp>

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace mixers {

template <class Type>
class broyden : public base<Type> {
	
public:

	broyden(const int arg_steps, const double arg_mix_factor, const long long dim, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
		iter_(0),
		max_size_(arg_steps),
		mix_factor_(arg_mix_factor),
		dv_({max_size_, dim}, NAN),
		df_({max_size_, dim}, NAN),
		f_old_(dim, NAN),
		vin_old_(dim, NAN),
		last_pos_(-1),
		comm_(comm){
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
		
	void broyden_extrapolation(math::array<Type, 1> & input_value, int const iter_used, math::array<Type, 1> const & ff){

		CALI_CXX_MARK_SCOPE("broyden_extrapolation");
		
		double const w0 = 0.01;
		double const ww = 5.0;

		if(iter_used == 0){
			gpu::run(input_value.size(),
							 [iv = begin(input_value), ffp = begin(ff),  mix = mix_factor_] GPU_LAMBDA (auto ip){
								 iv[ip] += mix*ffp[ip];
							 });

			return;
		}

		math::array<Type, 2> beta({iter_used, iter_used}, NAN);
		math::array<Type, 1> work(iter_used, NAN);

		//OPTIMIZATION: this should be done by gemm/gemv
		gpu::run(iter_used,
						 [iter_used, ivsize = input_value.size(), w0, ww, be = begin(beta), df = begin(df_), ffp = begin(ff), wo = begin(work), dfactor = 1.0/comm_.size()] GPU_LAMBDA (auto ii){
							 for(int jj = ii + 1; jj < iter_used; jj++){
								 Type aa = 0.0;
								 for(unsigned kk = 0; kk < ivsize; kk++) aa +=  ww*ww*conj(df[ii][kk])*df[jj][kk];
								 be[ii][jj] = aa;
								 be[jj][ii] = conj(aa);
							 }
							 be[ii][ii] = dfactor*(w0*w0 + ww*ww);

							 Type aa = 0.0;
							 for(unsigned kk = 0; kk < ivsize; kk++) aa += conj(df[ii][kk])*ffp[kk];
							 wo[ii] = aa;
						 });

		if(comm_.size() > 1){
			CALI_CXX_MARK_SCOPE("broyden_extrapolation::reduce");
			comm_.all_reduce_in_place_n(raw_pointer_cast(beta.data_elements()), beta.num_elements(), std::plus<>{});
			comm_.all_reduce_in_place_n(raw_pointer_cast(work.data_elements()), work.num_elements(), std::plus<>{});		
		}
		
		solvers::least_squares(beta, work);
		
		gpu::run(input_value.size(), 
						 [iv = begin(input_value),  ffp = begin(ff), ww, wo = begin(work), mix = mix_factor_, df = begin(df_), dv = begin(dv_), iter_used] GPU_LAMBDA (auto ip){
							 iv[ip] += mix*ffp[ip];							 
							 for(int ii = 0; ii < iter_used; ii++){
								 iv[ip] -= ww*ww*wo[ii]*(mix*df[ii][ip] + dv[ii][ip]);
							 }
						 });
			
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
		
	void operator()(math::array<Type, 1> & input_value, math::array<Type, 1> const & output_value){

		CALI_CXX_MARK_SCOPE("broyden_mixing");
		
		assert((typename math::array<double, 2>::size_type) input_value.size() == dv_[0].size());
		assert((typename math::array<double, 2>::size_type) output_value.size() == dv_[0].size());

		iter_++;

		math::array<Type, 1> ff(input_value.size());

		gpu::run(input_value.size(),
						 [iv = begin(input_value), ov = begin(output_value), ffp = begin(ff)] GPU_LAMBDA (auto ip){
								 ffp[ip] = ov[ip] - iv[ip]; 
						 });
		
		if(iter_ > 1){

			auto pos = (last_pos_ + 1)%max_size_;
				
			df_[pos] = ff;
			dv_[pos] = input_value;

			gpu::run(input_value.size(),
							 [pos, df = begin(df_), f_old = begin(f_old_), dv = begin(dv_), vin_old = begin(vin_old_)] GPU_LAMBDA (auto ip){
								 df[pos][ip] -= f_old[ip];
								 dv[pos][ip] -= vin_old[ip];
							 });

			using boost::multi::blas::dot;
			using boost::multi::blas::conj;
			
			gamma_ = dot(conj(df_[pos]), df_[pos]);

			if(comm_.size() > 1){
				CALI_CXX_MARK_SCOPE("broyden_mixing::reduce");
				comm_.all_reduce_in_place_n(&gamma_, 1, std::plus<>{});
			}

			gamma_ = std::max(1e-8, sqrt(gamma_));

			gpu::run(input_value.size(),
							 [pos, df = begin(df_), dv = begin(dv_), gamma = gamma_] GPU_LAMBDA (auto ip){
								 df[pos][ip] /= gamma;
								 dv[pos][ip] /= gamma;
							 });

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
	mutable boost::mpi3::communicator comm_;
	
};
	
}
}

#ifdef INQ_MIXERS_BROYDEN_UNIT_TEST
#undef INQ_MIXERS_BROYDEN_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("mixers::broyden", "[mixers::broyden]") {

	using namespace inq;
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
