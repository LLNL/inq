/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__ZERO
#define INQ__MATH__ZERO

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

#include <math/vector3.hpp>

#include <cassert>
#include <array>
#include <cmath>

#include <tinyformat/tinyformat.h>

namespace inq {

template <typename Type>
GPU_FUNCTION auto zero(){
  return Type{0.0};
}

template <>
GPU_FUNCTION auto zero<vector3<double>>(){
  return vector3<double>{0.0, 0.0, 0.0};
}

template <>
GPU_FUNCTION auto zero<vector3<double, covariant>>(){
  return vector3<double, covariant>{0.0, 0.0, 0.0};
}

template <>
GPU_FUNCTION auto zero<vector3<double, contravariant>>(){
  return vector3<double, contravariant>{0.0, 0.0, 0.0};
}

template <>
GPU_FUNCTION auto zero<vector3<complex>>(){
  return vector3<complex>{complex{0.0, 0.0}, complex{0.0, 0.0}, complex{0.0, 0.0}};
}

template <>
GPU_FUNCTION auto zero<vector3<complex, covariant>>(){
  return vector3<complex, covariant>{complex{0.0, 0.0}, complex{0.0, 0.0}, complex{0.0, 0.0}};
}

template <>
GPU_FUNCTION auto zero<vector3<complex, contravariant>>(){
  return vector3<complex, contravariant>{complex{0.0, 0.0}, complex{0.0, 0.0}, complex{0.0, 0.0}};
}

}
#endif

#ifdef INQ_MATH_ZERO_UNIT_TEST
#undef INQ_MATH_ZERO_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
}
#endif
