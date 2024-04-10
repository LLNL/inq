/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SCALAR_POTENTIAL
#define OPERATIONS__SCALAR_POTENTIAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/run.hpp>
#include <basis/field.hpp>
#include <states/orbital_set.hpp>

#include <cassert>

#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

template <class PotentialType, class ShiftType>
void scalar_potential_add(basis::field_set<basis::real_space, PotentialType> const & potential, int const index, ShiftType shift, states::orbital_set<basis::real_space, complex> const & phi, states::orbital_set<basis::real_space, complex> & vphi) {

	CALI_CXX_MARK_FUNCTION;
  
  assert(potential.basis() == phi.basis());
  
  gpu::run(phi.local_set_size(), phi.basis().local_size(),
           [pot = begin(potential.matrix()), it_vphi = begin(vphi.matrix()), it_phi = begin(phi.matrix()), shift, index] GPU_LAMBDA (auto ist, auto ip){
             it_vphi[ip][ist] += (pot[ip][index] + shift)*it_phi[ip][ist];
           });
}

template <class PotentialType, class ShiftType>
states::orbital_set<basis::real_space, complex> scalar_potential(basis::field_set<basis::real_space, PotentialType> const & potential, int const index, ShiftType shift, states::orbital_set<basis::real_space, complex> const & phi) {

	CALI_CXX_MARK_FUNCTION;

  states::orbital_set<basis::real_space, complex> vphi(phi.skeleton());
  
  assert(potential.basis() == phi.basis());
	
  gpu::run(phi.local_set_size(), phi.basis().local_size(),
           [pot = begin(potential.matrix()), it_vphi = begin(vphi.matrix()), it_phi = begin(phi.matrix()), shift, index] GPU_LAMBDA (auto ist, auto ip){
             it_vphi[ip][ist] = (pot[ip][index] + shift)*it_phi[ip][ist];
           });

  return vphi;
  
}

}
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_HAMILTONIAN_SCALAR_POTENTIAL_UNIT_TEST
#undef INQ_HAMILTONIAN_SCALAR_POTENTIAL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif

