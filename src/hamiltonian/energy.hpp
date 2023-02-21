/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__ENERGY
#define INQ__HAMILTONIAN__ENERGY

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

#include <tinyformat/tinyformat.h>

namespace inq {
namespace hamiltonian {

	struct energy {
		
		double ion;
		double ion_sr_lr;
		double eigenvalues;
		double external;
		double nonlocal;
		double hartree;
		double xc;
		double nvxc;
		double hf_exchange;

		energy(){
			ion = 0.0;
			eigenvalues = 0.0;
			external = 0.0;
			nonlocal = 0.0;
			hartree = 0.0;
			xc = 0.0;
			nvxc = 0.0;
			hf_exchange = 0.0;
		}

		template <typename HamType, typename ElType>
		auto calculate(HamType const & ham, ElType & el) {

			CALI_CXX_MARK_SCOPE("energy::calculate");

			auto normres = math::array<complex, 2>({el.lot().size(), el.max_local_size()});
			
			eigenvalues = 0.0;
			nonlocal = 0.0;
			hf_exchange = 0.0;
			
			int iphi = 0;
			for(auto & phi : el.lot()){
				
				auto residual = ham(phi);
				el.eigenvalues()[iphi] = operations::overlap_diagonal_normalized(residual, phi, operations::real_part{});
				operations::shift(-1.0, el.eigenvalues()[iphi], phi, residual);
				
				normres[iphi] = operations::overlap_diagonal(residual);
				auto nl_me = operations::overlap_diagonal_normalized(ham.non_local(phi), phi);
				auto exchange_me = operations::overlap_diagonal_normalized(ham.exchange(phi), phi);
				
				auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
				
				eigenvalues += operations::sum(el.occupations()[iphi], el.eigenvalues()[iphi], energy_term);
				nonlocal += operations::sum(el.occupations()[iphi], nl_me, energy_term);
				hf_exchange += 0.5*operations::sum(el.occupations()[iphi], exchange_me, energy_term);
				
				iphi++;
			}

			el.lot_states_comm_.all_reduce_n(&eigenvalues, 1);
			el.lot_states_comm_.all_reduce_n(&nonlocal, 1);
			el.lot_states_comm_.all_reduce_n(&hf_exchange, 1);

			return normres;
		}
	
		auto kinetic() const {
			return eigenvalues - 2.0*hartree - nvxc - 2.0*hf_exchange - external - nonlocal;
		}
		
		auto total() const {
			return kinetic() + hartree + external + nonlocal + xc + hf_exchange + ion;
		}

		template <class out_type>
		void print(out_type & out) const {

			tfm::format(out, "\n");
			tfm::format(out, "  total          = %20.12f\n", total());			
			tfm::format(out, "  kinetic        = %20.12f\n", kinetic());
			tfm::format(out, "  eigenvalues    = %20.12f\n", eigenvalues);
			tfm::format(out, "  hartree        = %20.12f\n", hartree);
			tfm::format(out, "  external       = %20.12f\n", external);
			tfm::format(out, "  nonlocal       = %20.12f\n", nonlocal);
			tfm::format(out, "  xc             = %20.12f\n", xc);
			tfm::format(out, "  intnvxc        = %20.12f\n", nvxc);
			tfm::format(out, "  HF exchange    = %20.12f\n", hf_exchange);
			tfm::format(out, "  ion            = %20.12f\n", ion);
			tfm::format(out, "\n");

		}
		
		template<class OStream>
		friend OStream& operator<<(OStream& os, energy const& self){
			self.print(os);
			return os;
		}
		
	};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ENERGY_UNIT_TEST
#undef INQ_HAMILTONIAN_ENERGY_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif

