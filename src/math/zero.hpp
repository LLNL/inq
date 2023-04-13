/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__ZERO
#define INQ__MATH__ZERO

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
