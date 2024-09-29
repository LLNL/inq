/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__FORCES
#define INQ__OBSERVABLES__FORCES

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <observables/density.hpp>
#include <operations/gradient.hpp>
#include <solvers/poisson.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace observables {

template <typename LongRangeType, typename ShortRangeType, typename GDensityType>
struct loc_pot {
	
	LongRangeType v1;
	ShortRangeType v2;
	GDensityType gdensityp;
	
	GPU_FUNCTION auto operator()(long ip) const {
		return (v1[ip] + v2[ip])*gdensityp[ip];
	}
};

struct forces_stress {
	gpu::array<vector3<double>, 1> forces;
	gpu::array<double, 2>          stress;

	forces_stress() = default;

	template <typename HamiltonianType>
	forces_stress(systems::ions const & ions, systems::electrons const & electrons, HamiltonianType const & ham):
		stress({3, 3}, 0.0)
	{
		forces = calculate_forces(ions, electrons, ham);
	}

private:
	
	template <typename HamiltonianType>
	gpu::array<vector3<double>, 1> calculate_forces(const systems::ions & ions, systems::electrons const & electrons, HamiltonianType const & ham){
		
		CALI_CXX_MARK_FUNCTION;
		
		basis::field<basis::real_space, vector3<double, covariant>> gdensity(electrons.density_basis());
		gdensity.fill(vector3<double, covariant>{0.0, 0.0, 0.0});
		
		gpu::array<vector3<double>, 1> forces_non_local(ions.size(), {0.0, 0.0, 0.0});
		
		auto iphi = 0;
		for(auto & phi : electrons.kpin()){
			
			auto gphi = operations::gradient(phi, /* factor = */ 1.0, /*shift = */ phi.kpoint());
			observables::density::calculate_gradient_add(electrons.occupations()[iphi], phi, gphi, gdensity);
			
			ham.projectors_all().force(phi, gphi, ions.cell().metric(), electrons.occupations()[iphi], ham.uniform_vector_potential(), forces_non_local);
			
			iphi++;
		}
		
		gdensity.all_reduce(electrons.kpin_states_comm());
		
		if(electrons.full_comm().size() > 1){
			CALI_CXX_MARK_SCOPE("forces_nonlocal::reduce");
			electrons.full_comm().all_reduce_n(raw_pointer_cast(forces_non_local.data_elements()), forces_non_local.size(), std::plus<>{});
		}
		
		auto ionic_forces = ionic::interaction_forces(ions.cell(), ions, electrons.atomic_pot());
		
		gpu::array<vector3<double>, 1> forces_local(ions.size(), {0.0, 0.0, 0.0});
		
		{ CALI_CXX_MARK_SCOPE("forces_local");
			
			solvers::poisson poisson_solver;
			
			//the force from the local potential
			for(int iatom = 0; iatom < ions.size(); iatom++){
				auto ionic_long_range = poisson_solver(electrons.atomic_pot().ionic_density(electrons.states_comm(), electrons.density_basis(), ions, iatom));
				auto ionic_short_range = electrons.atomic_pot().local_potential(electrons.states_comm(), electrons.density_basis(), ions, iatom);
				
				auto force_cov = -gpu::run(gpu::reduce(electrons.density_basis().local_size()),
																	 loc_pot<decltype(begin(ionic_long_range.linear())), decltype(begin(ionic_short_range.linear())), decltype(begin(gdensity.linear()))>
																	 {begin(ionic_long_range.linear()), begin(ionic_short_range.linear()), begin(gdensity.linear())});
				
				forces_local[iatom] = electrons.density_basis().volume_element()*ions.cell().metric().to_cartesian(force_cov);
			}
			
			if(electrons.density_basis().comm().size() > 1){
				CALI_CXX_MARK_SCOPE("forces_local::reduce");
				electrons.density_basis().comm().all_reduce_n(reinterpret_cast<double *>(raw_pointer_cast(forces_local.data_elements())), 3*forces_local.size());
			}
		}
		
		gpu::array<vector3<double>, 1> forces(ions.size());
		for(int iatom = 0; iatom < ions.size(); iatom++){
			forces[iatom] = ionic_forces[iatom] + forces_local[iatom] + forces_non_local[iatom];
		}
		
		// MISSING: the non-linear core correction term
		
		return forces;
	}
	
};

}
}
#endif

#ifdef INQ_OBSERVABLES_FORCES_UNIT_TEST
#undef INQ_OBSERVABLES_FORCES_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif

