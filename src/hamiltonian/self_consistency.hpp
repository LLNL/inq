/* -*- indent-tabs-mode: t -*- */

#ifndef HAMILTONIAN__KS_POTENTIAL
#define HAMILTONIAN__KS_POTENTIAL

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

#include <basis/field.hpp>
#include <functionals/lda.hpp>
#include <solvers/poisson.hpp>
#include <operations/sum.hpp>
#include <operations/integral.hpp>
#include <input/electronic_theory.hpp>

namespace hamiltonian {

	class self_consistency {

	public:

		self_consistency(input::electronic_theory arg_theory):
			theory_(arg_theory)	{

		}
		
		template <class vexternal_type, class density_type, class energy_type>
		auto ks_potential(const vexternal_type & vexternal, const density_type & electronic_density, const density_type & ionic_density, energy_type & energy){

			assert(vexternal.basis() == electronic_density.basis()); //for the moment they must be equal

			energy.external = operations::integral_product(electronic_density, vexternal);

			vexternal_type vks(vexternal.basis());

			solvers::poisson<basis::real_space> poisson_solver;

			switch(theory_){
				
			case input::electronic_theory::DENSITY_FUNCTIONAL:
				{

					auto total_density = operations::sum(electronic_density, ionic_density);
					auto vhartree = poisson_solver(total_density);
					energy.hartree = 0.5*operations::integral_product(total_density, vhartree);
					energy.nvhartree = operations::integral_product(electronic_density, vhartree);
						
					vexternal_type edxc(vexternal.basis());
					vexternal_type vxc(vexternal.basis());
					
					functionals::lda::xc_unpolarized(electronic_density.basis().size(), electronic_density, edxc, vxc);
					
					energy.xc = operations::integral_product(electronic_density, edxc);
					energy.nvxc = operations::integral_product(electronic_density, vxc);
					
					vks = operations::sum(vexternal, vhartree, vxc);

					break;
				}

			case input::electronic_theory::NON_INTERACTING:
				{

					auto vion = poisson_solver(ionic_density);
					energy.hartree = 0.5*operations::integral_product(ionic_density, vion);
					energy.nvhartree = operations::integral_product(electronic_density, vion);
					vks = operations::sum(vexternal, vion);
					
					break;
				}
				
			}
			
			return vks;			
		}

	private:

		input::electronic_theory theory_;

	};
}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::self_consistency", "[self_consistency]"){

  using namespace Catch::literals;
	
}

#endif

#endif
