/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SCALAR_POTENTIAL
#define OPERATIONS__SCALAR_POTENTIAL

/*
 Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa.

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

#include <gpu/run.hpp>
#include <basis/field.hpp>
#include <states/orbital_set.hpp>

#include <cassert>

#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

template <class PotentialType, class ShiftType>
void scalar_potential_add(basis::field<basis::real_space, PotentialType> const & potential, ShiftType shift, states::orbital_set<basis::real_space, complex> const & phi, states::orbital_set<basis::real_space, complex> & vphi) {

	CALI_CXX_MARK_FUNCTION;
  
  assert(potential.linear().num_elements() == phi.basis().local_size());
  
  gpu::run(phi.local_set_size(), phi.basis().local_size(),
           [pot = begin(potential.linear()), it_vphi = begin(vphi.matrix()), it_phi = begin(phi.matrix()), shift] GPU_LAMBDA (auto ist, auto ip){
             it_vphi[ip][ist] += (pot[ip] + shift)*it_phi[ip][ist];
           });
}

template <class PotentialType, class ShiftType>
states::orbital_set<basis::real_space, complex> scalar_potential(basis::field<basis::real_space, PotentialType> const & potential, ShiftType shift, states::orbital_set<basis::real_space, complex> const & phi) {

	CALI_CXX_MARK_FUNCTION;

  states::orbital_set<basis::real_space, complex> vphi(phi.skeleton());
  
  assert(potential.linear().num_elements() == phi.basis().local_size());
  
  gpu::run(phi.local_set_size(), phi.basis().local_size(),
           [pot = begin(potential.linear()), it_vphi = begin(vphi.matrix()), it_phi = begin(phi.matrix()), shift] GPU_LAMBDA (auto ist, auto ip){
             it_vphi[ip][ist] = (pot[ip] + shift)*it_phi[ip][ist];
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

