/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__DOUBLE_GRID
#define INQ__BASIS__DOUBLE_GRID

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>

#include <mpi3/environment.hpp>

#include <gpu/array.hpp>
#include <math/vector3.hpp>
#include <utils/interpolation_coefficients.hpp>
#include <parallel/partition.hpp>

namespace inq {
namespace basis {

class double_grid {
	
public:

  double_grid(bool enabled, int order = 5):
    min_(-order + 1),
    max_(order),
    enabled_(enabled)
  {

		if(not enabled_){
			min_ = 0;
			max_ = 0;
		}
		
    auto npoints = max_ - min_ + 1;
    gpu::array<double, 1> points(npoints);
    for(auto ii = 0; ii < npoints; ii++) points[ii] = min_ + ii;
    coeff_ = utils::interpolation_coefficients(points, 1.0/3.0);

  }

	struct reference {

	public:

		reference(int min, int max, gpu::array<double, 1>::const_iterator && coeff):
			min_(min),
			max_(max),
			coeff_(std::move(coeff))
		{
		}
		
		template <class Function>
		GPU_FUNCTION auto value(Function const & func, vector3<double> spacing, vector3<double> pos) const {
			
			if(min_ == 0) return func(pos);
			
			decltype(func(pos)) val = 0.0;
			
			for(int k0 = min_; k0 <= max_; k0++){
				for(int k1 = min_; k1 <= max_; k1++){
					for(int k2 = min_; k2 <= max_; k2++){
						
						double fac = coeff_[k0 - min_]*coeff_[k1 - min_]*coeff_[k2 - min_];
						
						for(int i0 = -1; i0 <= 1; i0++){
							for(int i1 = -1; i1 <= 1; i1++){
								for(int i2 = -1; i2 <= 1; i2++){
									val += fac*func(pos + spacing*vector3<double>{i0*(1.0/3.0 - k0), i1*(1.0/3.0 - k1), i2*(1.0/3.0 - k2)});
								}
							}
						}
						
						
					}
				}
			}
			
			return (1.0/27.0)*val;
		}

	private:
		
		int min_;
		int max_;
		gpu::array<double, 1>::const_iterator coeff_;
		
  };
		
	auto ref() const {
		return reference{min_, max_, begin(coeff_)};
	}
		
  auto spacing_factor() const {
    if(not enabled_) return 1.0;
    return 3.0;
  }

  auto & enabled() const {
    return enabled_;
  }
  
private:

  int min_;
  int max_;
  gpu::array<double, 1> coeff_;
  bool enabled_;
  
};


}
}

#endif

#ifdef INQ_BASIS_DOUBLE_GRID_UNIT_TEST
#undef INQ_BASIS_DOUBLE_GRID_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
	
	SECTION("disabled"){
		
		basis::double_grid dg(false);

		CHECK(not dg.enabled());
		CHECK(dg.spacing_factor() == 1.0);
		CHECK(dg.ref().value([](auto point){ return 1.0; }, {0.3, 0.3, 0.3}, {1.0, 2.0, 3.0})== 1.0_a);
		CHECK(dg.ref().value([](auto point){ return sin(point[1]); }, {0.3, 0.3, 0.3}, {1.0, 2.0, 3.0}) == 0.9092974268_a);

	}

	SECTION("enabled"){				
	
		basis::double_grid dg(true);
		
		CHECK(dg.enabled());
		CHECK(dg.spacing_factor() == 3.0);		
		CHECK(dg.ref().value([](auto point){ return 1.0; }, {0.3, 0.3, 0.3}, {1.0, 2.0, 3.0}) == 1.0_a);
		CHECK(dg.ref().value([](auto point){ return sin(point[1]); }, {0.3, 0.3, 0.3}, {1.0, 2.0, 3.0}) == 0.9092974261_a);
	}
  
}
#endif

