/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__KS_POTENTIAL
#define INQ__HAMILTONIAN__KS_POTENTIAL

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
#include <hamiltonian/atomic_potential.hpp>

#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

	class self_consistency {

	public:

		self_consistency(input::interaction interaction, basis::real_space const & potential_basis, basis::real_space const & density_basis):
			interaction_(interaction),
			exchange_(int(interaction.exchange())),
			correlation_(int(interaction.correlation())),
			vion_(density_basis),
			core_density_(density_basis),
			potential_basis_(potential_basis),
			density_basis_(density_basis)
		{
		}

		template <class ions_type>
		void update_ionic_fields(const ions_type & ions, const hamiltonian::atomic_potential & atomic_pot){

			CALI_CXX_MARK_FUNCTION;
			
			solvers::poisson poisson_solver;
			
			auto ionic_long_range = poisson_solver(atomic_pot.ionic_density(density_basis_, ions.cell(), ions.geo()));
			auto ionic_short_range = atomic_pot.local_potential(density_basis_, ions.cell(), ions.geo());
			vion_ = operations::add(ionic_long_range, ionic_short_range);

			core_density_ = atomic_pot.nlcc_density(density_basis_, ions.cell(), ions.geo());
		}
		
		template <class field_type, class energy_type>
		auto ks_potential(const field_type & electronic_density, energy_type & energy) const {

			CALI_CXX_MARK_FUNCTION;

			assert(electronic_density.basis() == density_basis_);
			assert(core_density_.basis() == electronic_density.basis());
			
			energy.external = operations::integral_product(electronic_density, vion_);

			field_type vks(vion_.skeleton());

			solvers::poisson poisson_solver;

			//IONIC POTENTIAL
			vks = vion_;

			// Hartree
			if(interaction_.hartree_potential()){
				auto vhartree = poisson_solver(electronic_density);
				energy.hartree = 0.5*operations::integral_product(electronic_density, vhartree);
				operations::increment(vks, vhartree);
			} else {
				energy.hartree = 0.0;
			}

			// XC
			energy.xc = 0.0;
			energy.nvxc = 0.0;
				
			if(exchange_.true_functional() or correlation_.true_functional()){

				auto full_density = operations::add(electronic_density, core_density_);
				double efunc = 0.0;
				field_type vfunc(vion_.skeleton());

				if(exchange_.true_functional()){
					exchange_(full_density, efunc, vfunc);
					energy.xc += efunc;
					operations::increment(vks, vfunc);
					energy.nvxc += operations::integral_product(electronic_density, vfunc); //the core correction does not go here
				}
				
				if(correlation_.true_functional()){
					correlation_(full_density, efunc, vfunc);
					energy.xc += efunc;
					operations::increment(vks, vfunc);
					energy.nvxc += operations::integral_product(electronic_density, vfunc); //the core correction does not go here
				}
			}
			
			if(potential_basis_ == vks.basis()){
				return vks;
			} else {
				return operations::transfer::coarsen(std::move(vks), potential_basis_);
			}
			
		}

	private:

		input::interaction interaction_;
		hamiltonian::xc_functional exchange_;
		hamiltonian::xc_functional correlation_;
		basis::field<basis::real_space, double> vion_;
		basis::field<basis::real_space, double> core_density_;
		basis::real_space potential_basis_;
		basis::real_space density_basis_;
		
	};
}
}

#ifdef INQ_HAMILTONIAN_SELF_CONSISTENCY_UNIT_TEST
#undef INQ_HAMILTONIAN_SELF_CONSISTENCY_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::self_consistency", "[self_consistency]"){

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif
