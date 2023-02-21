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
	calculate_energy(HamType const & ham, systems::electrons & el) :
		normres_    ({static_cast<size_type>(el.lot().size()), el.max_local_size()})
	{
			
		sum_eigenvalues_ = 0.0;
		nonlocal_ = 0.0;
		hf_exchange_ = 0.0;

		int iphi = 0;
		for(auto & phi : el.lot()){
			
			auto residual = ham(phi);
			el.eigenvalues()[iphi] = operations::overlap_diagonal_normalized(residual, phi, [] GPU_FUNCTION (auto xx){ return real(xx);});
			operations::shift(-1.0, el.eigenvalues()[iphi], phi, residual);
				
			normres_[iphi] = operations::overlap_diagonal(residual);
			auto nl_me = operations::overlap_diagonal_normalized(ham.non_local(phi), phi);
			auto exchange_me = operations::overlap_diagonal_normalized(ham.exchange(phi), phi);
				
			auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
				
			sum_eigenvalues_ += operations::sum(el.occupations()[iphi], el.eigenvalues()[iphi], energy_term);
			nonlocal_ += operations::sum(el.occupations()[iphi], nl_me, energy_term);
			hf_exchange_ += 0.5*operations::sum(el.occupations()[iphi], exchange_me, energy_term);

			iphi++;
		}

		el.lot_states_comm_.all_reduce_n(&sum_eigenvalues_, 1);
		el.lot_states_comm_.all_reduce_n(&nonlocal_       , 1);
		el.lot_states_comm_.all_reduce_n(&hf_exchange_    , 1);
	}
	
	math::array<complex, 2> normres_;

	using size_type = math::array<complex, 2>::size_type;
	using extensions_type = math::array<complex, 2>::extensions_type;

	double sum_eigenvalues_;
	double nonlocal_;
	double hf_exchange_;

};

}
}
#endif

#ifdef INQ_HAMILTONIAN_CALCULATE_ENERGY_UNIT_TEST
#undef INQ_HAMILTONIAN_CALCULATE_ENERGY_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif

