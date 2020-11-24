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

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <operations/gradient.hpp>
#include <solvers/poisson.hpp>

namespace inq {
namespace hamiltonian {

 math::array<math::vector3<double>, 1> calculate_forces(const systems::ions & ions, systems::electrons & electrons){
  solvers::poisson poisson_solver;
  
  math::array<math::vector3<double>, 1> forces(ions.geo().num_atoms(), {0.0, 0.0, 0.0});

  basis::field<basis::real_space, math::vector3<double>> gdensity(electrons.phi_.basis());

  //calculate the gradient of the density from the gradient of the orbitals  
  auto gphi = operations::gradient(electrons.phi_);
  
  gpu::run(3, electrons.phi_.basis().part().local_size(),
					 [nst = electrons.phi_.set_part().local_size(), occs = begin(electrons.states_.occupations()),
            phip = begin(electrons.phi_.matrix()), gphip = begin(gphi.matrix()), gdensityp = begin(gdensity.linear())] GPU_LAMBDA (auto idir, auto ip){
						 gdensityp[ip][idir] = 0.0;
						 for(int ist = 0; ist < nst; ist++) gdensityp[ip][idir] += occs[ist]*real(conj(gphip[ip][ist][idir])*phip[ip][ist] + conj(phip[ip][ist])*gphip[ip][ist][idir]);
					 });
  
	if(gphi.set_part().parallel()){
    gphi.set_comm().all_reduce_in_place_n(reinterpret_cast<double *>(static_cast<math::vec3d *>(gdensity.linear().data())), 3*gdensity.linear().size(), std::plus<>{});
	}

  //the force from the local potential  
  for(int iatom = 0; iatom < ions.geo().num_atoms(); iatom++){
    auto ionic_long_range = poisson_solver(electrons.atomic_pot_.ionic_density(electrons.density_basis_, ions.cell(), ions.geo(), iatom));
    auto ionic_short_range = electrons.atomic_pot_.local_potential(electrons.density_basis_, ions.cell(), ions.geo(), iatom);

    forces[iatom] = gpu::run(gpu::reduce(electrons.density_basis_.local_size()),
                             [v1 = begin(ionic_long_range.linear()), v2 = begin(ionic_short_range.linear()), gdensityp = begin(gdensity.linear())] GPU_LAMBDA (auto ip){
      return (v1[ip] + v2[ip])*gdensityp[ip];
    });
  }

  //the non-local potential term
  
  //the non-linear core correction term
  
  return forces;
}

}
}

#ifdef INQ_HAMILTONIAN_FORCES_UNIT_TEST
#undef INQ_HAMILTONIAN_FORCES_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::forces", "[forces]"){

	using namespace inq;
	using namespace Catch::literals;
	
}

#endif

#endif

