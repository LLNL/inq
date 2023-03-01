/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__CURRENT
#define INQ__OBSERVABLES__CURRENT

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <operations/gradient.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <physics/constants.hpp>

namespace inq {
namespace observables {

template <typename HamiltonianType>
basis::field<basis::real_space, vector3<double, covariant>> current_density(const systems::ions & ions, systems::electrons const & electrons, HamiltonianType const & ham){

	basis::field<basis::real_space, vector3<double, covariant>> cdensity(electrons.density_basis_);
	cdensity.fill(vector3<double, covariant>{0.0, 0.0, 0.0});

	auto iphi = 0;
	for(auto & phi : electrons.lot()){
		
		auto gphi = operations::gradient(phi, /* shift = */ ham.uniform_vector_potential(), /* factor = */ -1.0);		

		ham.projectors_all().position_commutator(phi, gphi, phi.kpoint() + ham.uniform_vector_potential());
		
    gpu::run(phi.basis().part().local_size(),
             [nst = phi.set_part().local_size(), occ = begin(electrons.occupations()[iphi]),
              ph = begin(phi.matrix()), gph = begin(gphi.matrix()), cdens = begin(cdensity.linear())] GPU_LAMBDA (auto ip){
               for(int ist = 0; ist < nst; ist++) cdens[ip] += 0.5*occ[ist]*imag(conj(ph[ip][ist])*gph[ip][ist] - conj(gph[ip][ist])*ph[ip][ist]);
             });
		iphi++;
	}
  
	cdensity.all_reduce(electrons.lot_states_comm_);
	return cdensity;
}

template <typename HamiltonianType>
auto current(const systems::ions & ions, systems::electrons const & electrons, HamiltonianType const & ham){
  return operations::integral(current_density(ions, electrons, ham));
}

}
}
#endif

#ifdef INQ_OBSERVABLES_CURRENT_UNIT_TEST
#undef INQ_OBSERVABLES_CURRENT_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
}
#endif
