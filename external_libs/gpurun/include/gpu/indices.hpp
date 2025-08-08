/* -*- indent-tabs-mode: t -*- */

#ifndef GPURUN__GPU__INDICES
#define GPURUN__GPU__INDICES

// Copyright (C) 2019-2025 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/host.hpp>

namespace gpu {

template <typename Int1, typename Int2, typename Int3>
GPU_FUNCTION auto bidimensional_to_linear(Int1 ix, Int2 iy, Int3 sizex) {
  return ix + iy*sizex;
}

template <typename Int1, typename Int2, typename Int3, typename Int4>
GPU_FUNCTION void linear_to_bidimensional(Int1 ind, Int2 sizex, Int3 & ix, Int4 & iy) {
  ix = ind%sizex;
  iy = ind/sizex;
}

template <typename Int1, typename Int2, typename Int3, typename Int4, typename Int5>
GPU_FUNCTION auto tridimensional_to_linear(Int1 ix, Int2 iy, Int3 iz, Int4 sizex, Int5 sizey) {
  return ix + (iy + iz*sizey)*sizex;
}

template <typename Int1, typename Int2, typename Int3, typename Int4, typename Int5, typename Int6>
GPU_FUNCTION void linear_to_tridimensional(Int1 ind, Int2 sizex, Int3 sizey, Int4 & ix, Int5 & iy, Int6 & iz) {
  ix = ind%sizex;
  iy = (ind/sizex)%sizey;
  iz = ind/(sizex*sizey);
}

}
#endif

#ifdef GPURUN__INDICES__UNIT_TEST
#undef GPURUN__INDICES__UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(GPURUN_TEST_FILE, GPURUN_TEST_TAG) {

  using namespace gpu;
  
  SECTION("2D") {
  
    auto nii = 456;
    auto njj = 519;
    
    auto idx = 0;
    for(auto jj = 0; jj < njj; jj++) {
      for(auto ii = 0; ii < nii; ii++) {
        
        int ii2, jj2;
        
        CHECK(bidimensional_to_linear(ii, jj, nii) == idx);
        linear_to_bidimensional(idx, nii, ii2, jj2);
        CHECK(ii2 == ii);
        CHECK(jj2 == jj);
        
        idx++;
      }
    }
    CHECK(idx == nii*njj);
  }

    SECTION("2D") {
  
    auto nii = 456;
    auto njj = 519;
    auto nkk = 347;
    
    auto idx = 0;
    for(auto kk = 0; kk < nkk; kk++) {
      for(auto jj = 0; jj < njj; jj++) {
        for(auto ii = 0; ii < nii; ii++) {
          
          int ii2, jj2, kk2;
          
          CHECK(tridimensional_to_linear(ii, jj, kk, nii, njj) == idx);
          linear_to_tridimensional(idx, nii, njj, ii2, jj2, kk2);
          CHECK(ii2 == ii);
          CHECK(jj2 == jj);
          CHECK(kk2 == kk);          
          
          idx++;
        }
      }
    }
    
    CHECK(idx == nii*njj*nkk);
  }
  
  
}
#endif
