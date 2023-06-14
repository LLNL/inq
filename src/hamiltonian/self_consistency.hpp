/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__KS_POTENTIAL
#define INQ__HAMILTONIAN__KS_POTENTIAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <basis/field.hpp>
#include <solvers/poisson.hpp>
#include <observables/density.hpp>
#include <operations/add.hpp>
#include <operations/integral.hpp>
#include <input/interaction.hpp>
#include <hamiltonian/xc_term.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <perturbations/none.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

template <typename Perturbation = perturbations::none>
class self_consistency {
	
public:
	
	self_consistency(input::interaction interaction, basis::real_space const & potential_basis, basis::real_space const & density_basis, int const spin_components, Perturbation const & pert = {}):
		interaction_(interaction),
		xc_(interaction, spin_components),
		vion_(density_basis),
		core_density_(density_basis),
		potential_basis_(potential_basis),
		density_basis_(density_basis),
		pert_(pert)
	{
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////
	
	self_consistency(self_consistency && old, parallel::communicator new_comm):
		interaction_(std::move(old.interaction_)),
		xc_(std::move(old.xc_)),
		vion_(std::move(old.vion_), new_comm),
		core_density_(std::move(old.core_density_), new_comm),
		potential_basis_(std::move(old.potential_basis_), new_comm),
		density_basis_(std::move(old.density_basis_), new_comm),
		pert_(std::move(old.pert_))
	{
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////
	
	template <class CommType, class ions_type>
	void update_ionic_fields(CommType & comm, const ions_type & ions, const hamiltonian::atomic_potential & atomic_pot){
		
		CALI_CXX_MARK_FUNCTION;
		
		solvers::poisson poisson_solver;
		
		auto ionic_long_range = poisson_solver(atomic_pot.ionic_density(comm, density_basis_, ions.geo()));
		auto ionic_short_range = atomic_pot.local_potential(comm, density_basis_, ions.geo());
		vion_ = operations::add(ionic_long_range, ionic_short_range);
		
		core_density_ = atomic_pot.nlcc_density(comm, density_basis_, ions.geo());
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	
	template <typename HamiltonianType, typename EnergyType, typename FieldType>
	void update_hamiltonian(HamiltonianType & hamiltonian, EnergyType & energy, FieldType const & spin_density, double time = 0.0) const {

		CALI_CXX_MARK_FUNCTION;

		assert(spin_density.basis() == density_basis_);
		assert(core_density_.basis() == spin_density.basis());

		auto total_density = observables::density::total(spin_density);
			
		energy.external(operations::integral_product(total_density, vion_));

		std::vector<basis::field<basis::real_space, typename HamiltonianType::potential_type>> vks(spin_density.set_size(), density_basis_);

		solvers::poisson poisson_solver;

		//IONIC POTENTIAL
		for(auto & vcomp : vks) vcomp = vion_;

		//Time-dependent perturbation
		if(pert_.has_uniform_electric_field()){
			auto efield = pert_.uniform_electric_field(time);

			for(auto & vcomp : vks) {
				gpu::run(vcomp.basis().local_sizes()[2], vcomp.basis().local_sizes()[1], vcomp.basis().local_sizes()[0],
								 [point_op = vcomp.basis().point_op(), efield, vk = begin(vcomp.cubic())] GPU_LAMBDA (auto iz, auto iy, auto ix){
									 auto rr = point_op.rvector_cartesian(ix, iy, iz);
									 vk[ix][iy][iz] += -dot(efield, rr);
								 });
			}
		}

		for(auto & vcomp : vks) {
			pert_.potential(time, vcomp);
		}
		
		// Hartree
		if(interaction_.hartree_potential()){
			auto vhartree = poisson_solver(total_density);
			energy.hartree(0.5*operations::integral_product(total_density, vhartree));
			for(auto & vcomp : vks) operations::increment(vcomp, vhartree);
		} else {
			energy.hartree(0.0);
		}

		// XC
		double exc, nvxc;
		xc_(spin_density, core_density_, vks, exc, nvxc);
		energy.xc(exc);
		energy.nvxc(nvxc);

		// PUT THE CALCULATED POTENTIAL IN THE HAMILTONIAN
		if(potential_basis_ == vks[0].basis()){
			hamiltonian.scalar_potential_= std::move(vks);
		} else {
			for(unsigned ipot = 0; ipot < vks.size(); ipot++){
				hamiltonian.scalar_potential_[ipot] = operations::transfer::coarsen(vks[ipot], potential_basis_);
			}
		}

		// THE VECTOR POTENTIAL
		
		if(pert_.has_uniform_vector_potential()){
			hamiltonian.uniform_vector_potential_ = potential_basis_.cell().metric().to_covariant(pert_.uniform_vector_potential(time));
		} else {
			hamiltonian.uniform_vector_potential_ = {0.0, 0.0, 0.0};
		}
		
	}

	void update_induced_potential(vector3<double,covariant> & induced, vector3<double,covariant> & velocity, vector3<double,covariant> & accel, double const dt, const double volume, vector3<double,covariant> const & current) const {
		if(interaction_.induced_vector_potential_value() == input::interaction::induced_vector_potential::GAUGE_FIELD){
			induced += 0.5*dt*velocity;
			velocity += 0.5*dt*accel;
			accel = interaction_.alpha_value()*(-current)/volume;
			velocity += 0.5*dt*accel;
			induced += 0.5*dt*velocity;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	
	auto exx_coefficient(){
		if(xc_.exchange().true_functional()) return xc_.exchange().exx_coefficient();
		return interaction_.exchange_coefficient();
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	
private:
	
	input::interaction interaction_;
	hamiltonian::xc_term xc_;
	basis::field<basis::real_space, double> vion_;
	basis::field<basis::real_space, double> core_density_;
	basis::real_space potential_basis_;
	basis::real_space density_basis_;
	Perturbation pert_;
		
};
}
}
#endif

#ifdef INQ_HAMILTONIAN_SELF_CONSISTENCY_UNIT_TEST
#undef INQ_HAMILTONIAN_SELF_CONSISTENCY_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;

}
#endif
