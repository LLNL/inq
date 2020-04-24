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
#include <solvers/poisson.hpp>
#include <operations/add.hpp>
#include <operations/integral.hpp>
#include <input/interaction.hpp>
#include <hamiltonian/xc_functional.hpp>

namespace hamiltonian {

	class self_consistency {

	public:

		self_consistency(input::interaction interaction):
			theory_(interaction.theory()),
			exchange_(int(interaction.exchange())),
			correlation_(int(interaction.correlation())){
		}
		
		template <class field_type, class energy_type>
		auto ks_potential(const field_type & vexternal, const field_type & electronic_density, const field_type & core_density, const field_type & ionic_density, energy_type & energy){

			assert(vexternal.basis() == electronic_density.basis()); //for the moment they must be equal

			energy.external = operations::integral_product(electronic_density, vexternal);

			field_type vks(vexternal.skeleton());

			solvers::poisson poisson_solver;

			switch(theory_){

			case input::interaction::electronic_theory::HARTREE_FOCK:
				{

					auto total_density = operations::add(electronic_density, ionic_density);
					auto vhartree = poisson_solver(total_density);
					energy.hartree = 0.5*operations::integral_product(electronic_density, vhartree);
					
					vks = operations::add(vexternal, vhartree);
					
					break;
				}
				
			case input::interaction::electronic_theory::DENSITY_FUNCTIONAL:
				{

					auto vhartree = poisson_solver(electronic_density);
					auto vion = poisson_solver(ionic_density);
					
					energy.hartree = 0.5*operations::integral_product(electronic_density, vhartree);
					energy.external += operations::integral_product(electronic_density, vion);					
					vion = operations::add(vion, vexternal);
					
					double ex, ec;
					field_type vx(vexternal.skeleton());
					field_type vc(vexternal.skeleton());

					{
						auto full_density = operations::add(electronic_density, core_density);
						exchange_(full_density, ex, vx);
						correlation_(full_density, ec, vc);
					}

					energy.xc = ex + ec;
					auto vxc = operations::add(vx, vc);
					energy.nvxc = operations::integral_product(electronic_density, vxc); //the core correction does not go here

					vks = operations::add(vion, vhartree, vxc);

					break;
				}

			case input::interaction::electronic_theory::NON_INTERACTING:
				{

					auto vion = poisson_solver(ionic_density);
					energy.external += operations::integral_product(electronic_density, vion);
					vks = operations::add(vexternal, vion);
					
					break;
				}
				
			}
			
			return vks;
		}

		auto theory() const {
			return theory_;
		}

	private:

		input::interaction::electronic_theory theory_;
		hamiltonian::xc_functional exchange_;
		hamiltonian::xc_functional correlation_;

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
