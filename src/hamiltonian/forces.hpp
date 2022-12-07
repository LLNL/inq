/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__FORCES
#define INQ__HAMILTONIAN__FORCES

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

#include <observables/density.hpp>
#include <operations/gradient.hpp>
#include <solvers/poisson.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {

template <typename LongRangeType, typename ShortRangeType, typename GDensityType>
struct loc_pot {
	
	LongRangeType v1;
	ShortRangeType v2;
	GDensityType gdensityp;
	
	GPU_FUNCTION auto operator()(long ip) const {
		return (v1[ip] + v2[ip])*gdensityp[ip];
	}

};


template <typename HamiltonianType>
math::array<math::vector3<double>, 1> calculate_forces(const systems::ions & ions, systems::electrons & electrons, HamiltonianType const & ham){

	CALI_CXX_MARK_FUNCTION;

	basis::field<basis::real_space, math::vector3<double, math::covariant>> gdensity(electrons.density_basis_);
	gdensity = {0.0, 0.0, 0.0};

  math::array<math::vector3<double>, 1> forces_non_local(ions.geo().num_atoms(), {0.0, 0.0, 0.0});

	auto iphi = 0;
	for(auto & phi : electrons.lot()){
		
		auto gphi = operations::gradient(phi);
		observables::density::calculate_gradient_add(electrons.occupations()[iphi], phi, gphi, gdensity);

		ham.projectors_all().force(phi, gphi, ions.cell().metric(), electrons.occupations()[iphi], ham.uniform_vector_potential(), forces_non_local);
		
		iphi++;
	}
	
	if(electrons.lot_states_comm_.size() > 1){
		CALI_CXX_MARK_SCOPE("forces::gdensity::reduce");
		electrons.lot_states_comm_.all_reduce_in_place_n(reinterpret_cast<double *>(raw_pointer_cast(gdensity.linear().data_elements())), 3*gdensity.linear().size(), std::plus<>{});
	}
	
	if(electrons.full_comm_.size() > 1){
		CALI_CXX_MARK_SCOPE("forces_nonlocal::reduce");
		electrons.full_comm_.all_reduce_in_place_n(reinterpret_cast<double *>(raw_pointer_cast(forces_non_local.data_elements())), 3*forces_non_local.size(), std::plus<>{});
	}
	
	//ionic force
	auto ionic_forces = inq::ions::interaction_forces(ions.cell(), ions.geo(), electrons.atomic_pot_);

	math::array<math::vector3<double>, 1> forces_local(ions.geo().num_atoms(), {0.0, 0.0, 0.0});

	{
		CALI_CXX_MARK_SCOPE("forces_local");
		
		solvers::poisson poisson_solver;
		
		//the force from the local potential
		for(int iatom = 0; iatom < ions.geo().num_atoms(); iatom++){
			auto ionic_long_range = poisson_solver(electrons.atomic_pot_.ionic_density(electrons.states_comm_, electrons.density_basis_, ions.geo(), iatom));
			auto ionic_short_range = electrons.atomic_pot_.local_potential(electrons.states_comm_, electrons.density_basis_, ions.geo(), iatom);

			auto force_cov = -gpu::run(gpu::reduce(electrons.density_basis_.local_size()),
																 loc_pot<decltype(begin(ionic_long_range.linear())), decltype(begin(ionic_short_range.linear())), decltype(begin(gdensity.linear()))>
																 {begin(ionic_long_range.linear()), begin(ionic_short_range.linear()), begin(gdensity.linear())});
			
			forces_local[iatom] = electrons.density_basis_.volume_element()*ions.cell().metric().to_cartesian(force_cov);
		}

		if(electrons.density_basis_.comm().size() > 1){
			CALI_CXX_MARK_SCOPE("forces_local::reduce");
			electrons.density_basis_.comm().all_reduce_in_place_n(reinterpret_cast<double *>(raw_pointer_cast(forces_local.data_elements())), 3*forces_local.size(), std::plus<>{});
		}
		
	}

	math::array<math::vector3<double>, 1> forces(ions.geo().num_atoms());
  for(int iatom = 0; iatom < ions.geo().num_atoms(); iatom++){
    forces[iatom] = ionic_forces[iatom] + forces_local[iatom] + forces_non_local[iatom];
  }
	
  // MISSING: the non-linear core correction term
  
  return forces;
}

}
}

#ifdef INQ_HAMILTONIAN_FORCES_UNIT_TEST
#undef INQ_HAMILTONIAN_FORCES_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::forces", "[forces]"){

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif

