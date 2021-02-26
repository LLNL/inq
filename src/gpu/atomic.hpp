/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GPU__ATOMIC
#define INQ__GPU__ATOMIC

/*
 Copyright (C) 2021 Xavier Andrade

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

namespace inq {
namespace gpu {
namespace atomic {

template <typename Type1, typename Type2>
GPU_FUNCTION inline Type1 add(Type1 * val, Type2 const & incr){
#ifdef ENABLE_CUDA
  return atomicAdd(val, incr);
#else
  auto old_val = *val;
  *val += incr;
  return old_val;
#endif
}

#ifdef ENABLE_CUDA
template <typename Type2>
GPU_FUNCTION inline long add(long * val, Type2 const & incr){
  return add((unsigned long long int *) val, incr);
}
#endif

}
}
}

#ifdef INQ_GPU_ATOMIC_UNIT_TEST
#undef INQ_GPU_ATOMIC_UNIT_TEST

#endif
#endif
