/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MIXERS__BROYDEN
#define INQ__MIXERS__BROYDEN

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <math/complex.hpp>
#include <math/vector3.hpp>
#include <gpu/array.hpp>
#include <solvers/least_squares.hpp>
#include <mixers/base.hpp>

#include <parallel/communicator.hpp>
#include <mpi3/environment.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace mixers {

template <class ArrayType>
class broyden : public base<ArrayType> {
	
public:

	using element_type = typename ArrayType::element_type;
	
	template <class CommType>
	broyden(const int arg_steps, const double arg_mix_factor, const long long dim, CommType & comm):
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
		
	void broyden_extrapolation(ArrayType & input_value, int const iter_used, gpu::array<element_type, 1> const & ff){

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

		namespace blas = boost::multi::blas;
		
		auto subdf = df_({0, iter_used}, {0, input_value.size()});
		auto beta = +blas::gemm(ww*ww, subdf, blas::H(subdf));

		auto matff = gpu::array<element_type, 2>({iter_used, input_value.size()}); //for some reason this has to be iter_used and not 1. 
		matff[0] = ff({0, input_value.size()});
		auto workmat = +blas::gemm(1.0, matff, blas::H(subdf));
		auto work = +workmat[0];
		
		gpu::run(iter_used, [w0, ww, be = begin(beta), dfactor = 1.0/comm_.size()] GPU_LAMBDA (auto ii){ be[ii][ii] = dfactor*(w0*w0 + ww*ww); });
		
		if(comm_.size() > 1){
			CALI_CXX_MARK_SCOPE("broyden_extrapolation::reduce");
			comm_.all_reduce_n(raw_pointer_cast(beta.data_elements()), beta.num_elements());
			comm_.all_reduce_n(raw_pointer_cast(work.data_elements()), work.num_elements());
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
		
	void operator()(ArrayType & input_value, ArrayType const & output_value){

		CALI_CXX_MARK_SCOPE("broyden_mixing");
		
		assert((typename gpu::array<double, 2>::size_type) input_value.size() == dv_[0].size());
		assert((typename gpu::array<double, 2>::size_type) output_value.size() == dv_[0].size());

		iter_++;

		gpu::array<element_type, 1> ff(input_value.size());

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
	gpu::array<element_type, 2> dv_;
	gpu::array<element_type, 2> df_;
	gpu::array<element_type, 1> f_old_;
	gpu::array<element_type, 1> vin_old_;
	element_type gamma_;
	int last_pos_;
	mutable parallel::communicator comm_;
	
};
	
}
}
#endif

#ifdef INQ_MIXERS_BROYDEN_UNIT_TEST
#undef INQ_MIXERS_BROYDEN_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
 

	gpu::array<double, 1> vin({10.0, -20.0});
	gpu::array<double, 1> vout({0.0,  22.2});

  mixers::broyden<decltype(vin)> lm(5, 0.5, 2, boost::mpi3::environment::get_self_instance());
	
	lm(vin, vout);
  
	CHECK(vin[0] == 5.0_a);
  CHECK(vin[1] == 1.1_a);

	vout = gpu::array<double, 1>({4.0, 5.5});

	lm(vin, vout);

	CHECK(vin[0] == 4.4419411001_a);
  CHECK(vin[1] == 3.5554591594_a);

}
#endif
