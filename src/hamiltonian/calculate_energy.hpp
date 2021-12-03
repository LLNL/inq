/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__CALCULATE_ENERGY
#define INQ__HAMILTONIAN__CALCULATE_ENERGY

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

#include <operations/shift.hpp>
#include <operations/overlap_diagonal.hpp>
#include <systems/electrons.hpp>

namespace inq {
namespace hamiltonian {

  class calculate_energy{

  public:

    template <typename HamType>
    calculate_energy(HamType const & ham, systems::electrons const & el){
      
      auto residual = ham(el.phi());
      eigenvalues_ = operations::overlap_diagonal_normalized(residual, el.phi());
      operations::shift(-1.0, eigenvalues_, el.phi().fields(), residual);
      
      normres_ = operations::overlap_diagonal(residual);
      auto nl_me = operations::overlap_diagonal_normalized(ham.non_local(el.phi()), el.phi());
      auto exchange_me = operations::overlap_diagonal_normalized(ham.exchange(el.phi().fields()), el.phi().fields());
      
      auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
      
      sum_eigenvalues_ = operations::sum(el.occupations()[0], eigenvalues_, energy_term);
      nonlocal_ = operations::sum(el.occupations()[0], nl_me, energy_term);
      hf_exchange_ = operations::sum(el.occupations()[0], exchange_me, energy_term);

      state_conv_ = operations::sum(el.occupations()[0], normres_, [](auto occ, auto nres){ return fabs(occ)*fabs(nres); });
					
			state_conv_ /= el.states_.num_electrons();
      
      el.phi().fields().set_comm().all_reduce_in_place_n(&sum_eigenvalues_, 1, std::plus<>{});
      el.phi().fields().set_comm().all_reduce_in_place_n(&nonlocal_, 1, std::plus<>{});
      el.phi().fields().set_comm().all_reduce_in_place_n(&hf_exchange_, 1, std::plus<>{});
      el.phi().fields().set_comm().all_reduce_in_place_n(&state_conv_, 1, std::plus<>{});
      
    }

    math::array<complex, 1> normres_;
    math::array<complex, 1> eigenvalues_;
    double sum_eigenvalues_;
    double nonlocal_;
    double hf_exchange_;
    double state_conv_;
  };
  
}
}

#ifdef INQ_HAMILTONIAN_CALCULATE_ENERGY_UNIT_TEST
#undef INQ_HAMILTONIAN_CALCULATE_ENERGY_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::calculate_energy", "[calculate_energy]"){

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif

