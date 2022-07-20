/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GPU__COPY
#define INQ__GPU__COPY

/*
 Copyright (C) 2022 Xavier Andrade

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

#include <inq_config.h>

#include <gpu/run.hpp>
#include <math/complex.hpp>

namespace inq {
namespace gpu {

template <typename SourceType, typename DestinationType>
void copy(long dim1, long dim2, SourceType const & source, DestinationType & destination){
  gpu::run(dim2, dim1,
           [sou = begin(source), des = begin(destination)] GPU_LAMBDA (auto i2, auto i1){
             des[i1][i2] = sou[i1][i2];
           });
}

}
}

#ifdef INQ_GPU_COPY_UNIT_TEST
#undef INQ_GPU_COPY_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <math/array.hpp>

TEST_CASE("gpu::copy", "[gpu::copy]") {
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
#endif
