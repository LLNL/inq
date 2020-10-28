/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GPU__REDUCE
#define INQ__GPU__REDUCE

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

#include <inq_config.h>

#ifdef ENABLE_CUDA
#include <cuda.h>
#endif

#include <cstddef> // std::size_t
#include <cassert>

namespace inq {
namespace gpu {

template <class kernel_type>
auto reduce(size_t size, kernel_type kernel){

  using type = decltype(kernel(0));
  
#ifndef ENABLE_CUDA

  type accumulator = 0.0;
  for(size_t ii = 0; ii < size; ii++){
    accumulator += kernel(ii);
  }
  return accumulator;

#else

#endif

  
  
}

}
}

#ifdef INQ_GPU_REDUCE_UNIT_TEST
#undef INQ_GPU_REDUCE_UNIT_TEST

#include <catch2/catch.hpp>
#include <math/array.hpp>

#include <gpu/run.hpp>

TEST_CASE("function gpu::reduce", "[gpu::reduce]") {

	using namespace inq;
	using namespace Catch::literals;

  const size_t maxsize = 387420490;

  for(size_t nn = 1; nn < maxsize; nn *= 3){
    CHECK(gpu::reduce(nn, [] GPU_LAMBDA (auto ii){ return double(ii); }) == (nn*(nn - 1.0)/2.0));
  }
  
}

#endif
#endif
