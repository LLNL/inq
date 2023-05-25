/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATRIX__DIAGONAL
#define INQ__MATRIX__DIAGONAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <matrix/gather_scatter.hpp>

namespace inq {
namespace matrix {

template <typename DistributedMatrix>
auto diagonal(DistributedMatrix const & matrix) -> gpu::array<typename DistributedMatrix::element_type, 1> {

  assert(matrix.sizex() == matrix.sizey());

  auto part = matrix.party();
  auto array = matrix::all_gather(matrix);
  
  gpu::array<typename DistributedMatrix::element_type, 1> diag(part.local_size());
  
  gpu::run(part.local_size(), [dia = begin(diag), arr = begin(array), pa = part] GPU_LAMBDA (auto ii){
    auto iig = pa.start() + ii;
    dia[ii] = arr[iig][iig];
  });
  
  return diag;
}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_MATRIX_DIAGONAL_UNIT_TEST
#undef INQ_MATRIX_DIAGONAL_UNIT_TEST

#include <gpu/array.hpp>
#include <matrix/distributed.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

}

#endif
