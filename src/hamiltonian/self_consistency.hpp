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
#include <perturbations/none.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

template <typename Perturbation = perturbations::none>
	class self_consistency {

	public:

	self_consistency(input::interaction interaction, basis::real_space const & potential_basis, basis::real_space const & density_basis, Perturbation const & pert = {}):
			interaction_(interaction),
			exchange_(int(interaction.exchange())),
			correlation_(int(interaction.correlation())),
			vion_(density_basis),
			core_density_(density_basis),
			potential_basis_(potential_basis),
			density_basis_(density_basis),
			pert_(pert)
		{
		}

		self_consistency(self_consistency && old, parallel::communicator new_comm):
			interaction_(std::move(old.interaction_)),
			exchange_(int(interaction_.exchange())),
			correlation_(int(interaction_.correlation())),
			vion_(std::move(old.vion_), new_comm),
			core_density_(std::move(old.core_density_), new_comm),
			potential_basis_(std::move(old.potential_basis_), new_comm),
			density_basis_(std::move(old.density_basis_), new_comm),
			pert_(std::move(old.pert_))
		{
		}
				
		template <class CommType, class ions_type>
		void update_ionic_fields(CommType & comm, const ions_type & ions, const hamiltonian::atomic_potential & atomic_pot){

			CALI_CXX_MARK_FUNCTION;
			
			solvers::poisson poisson_solver;
			
			auto ionic_long_range = poisson_solver(atomic_pot.ionic_density(comm, density_basis_, ions.cell(), ions.geo()));
			auto ionic_short_range = atomic_pot.local_potential(comm, density_basis_, ions.cell(), ions.geo());
			vion_ = operations::add(ionic_long_range, ionic_short_range);

			core_density_ = atomic_pot.nlcc_density(comm, density_basis_, ions.cell(), ions.geo());
		}
		
		template <class field_type, class energy_type>
		field_type ks_potential(const field_type & electronic_density, energy_type & energy, double time = 0.0) const {

			CALI_CXX_MARK_FUNCTION;

			assert(electronic_density.basis() == density_basis_);
			assert(core_density_.basis() == electronic_density.basis());
			
			energy.external = operations::integral_product(electronic_density, vion_);

			field_type vks(vion_.skeleton());

			solvers::poisson poisson_solver;

			//IONIC POTENTIAL
			vks = vion_;

			//Time-dependent perturbation
			if(pert_.has_uniform_electric_field()){
				auto efield = pert_.uniform_electric_field(time);
				gpu::run(vks.basis().local_sizes()[2], vks.basis().local_sizes()[1], vks.basis().local_sizes()[0],
								 [point_op = vks.basis().point_op(), efield, vk = begin(vks.cubic())] GPU_LAMBDA (auto iz, auto iy, auto ix){
									 auto rr = point_op.rvector_cartesian(ix, iy, iz);
									 vk[ix][iy][iz] += -dot(efield, rr);
								 });
			}
			
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

		auto exx_coefficient(){
			if(exchange_.true_functional()) return exchange_.exx_coefficient();
			return interaction_.exchange_coefficient();
		}

	private:

		input::interaction interaction_;
		hamiltonian::xc_functional exchange_;
		hamiltonian::xc_functional correlation_;
		basis::field<basis::real_space, double> vion_;
		basis::field<basis::real_space, double> core_density_;
		basis::real_space potential_basis_;
		basis::real_space density_basis_;
	Perturbation pert_;
		
	};
}
}

#ifdef INQ_HAMILTONIAN_SELF_CONSISTENCY_UNIT_TEST
#undef INQ_HAMILTONIAN_SELF_CONSISTENCY_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::self_consistency", "[self_consistency]"){

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif
