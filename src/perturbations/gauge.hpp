/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__GAUGE
#define INQ__PERTURBATIONS__GAUGE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <states/orbital_set.hpp>

namespace inq {
namespace perturbations {

enum class gauge {
  mixed,    // length gauge in non-periodic dimensions, velocity in periodic dimensions
  length,   // the electric field enters through the scalar potential or a phase in the orbitals, does not work for periodic dimensions
  velocity  // the electric field is applied through a uniform vector potential
};

}
}
#endif

#ifdef INQ_PERTURBATIONS_GAUGE_UNIT_TEST
#undef INQ_PERTURBATIONS_GAUGE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  [[maybe_unused]] auto gau = perturbations::gauge::velocity;

}
#endif
