/* -*- indent-tabs-mode: t -*- */

#ifndef GPURUN__GPU__COPY
#define GPURUN__GPU__COPY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/run.hpp>
#include <math/complex.hpp>

namespace gpu {

template <typename SourceType, typename DestinationType>
void copy(long dim1, long dim2, SourceType const & source, DestinationType & destination){
  gpu::run(dim2, dim1,
           [sou = begin(source), des = begin(destination)] GPU_LAMBDA (auto i2, auto i1){
             des[i1][i2] = sou[i1][i2];
           });
}

}
#endif

#ifdef GPURUN__COPY__UNIT_TEST
#undef GPURUN__COPY__UNIT_TEST

#include <catch2/catch_all.hpp>
#include <math/array.hpp>

TEST_CASE(GPURUN_TEST_FILE, GPURUN_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
  
  SECTION("2D"){
    long nn = 1000;
    long mm = 450;
    
    math::array<double, 2> src({nn + 2, mm + 2}, 10.0);
    math::array<double, 2> dest({nn, mm}, 0.0);

    for(int in = 0; in < nn; in++){
      for(int im = 0; im < mm; im++){
        src[in][im] = in + nn*im;
      }
    }
    
    gpu::copy(nn, mm, src, dest);

    for(int in = 0; in < nn; in++){
      for(int im = 0; im < mm; im++){
        CHECK(dest[in][im] == in + nn*im);        
      }
    }
  }
  
}
#endif
