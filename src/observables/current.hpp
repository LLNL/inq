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
basis::field<basis::real_space, math::vector3<double, math::covariant>> current_density(const systems::ions & ions, systems::electrons & electrons, HamiltonianType const & ham){

	basis::field<basis::real_space, math::vector3<double, math::covariant>> cdensity(electrons.density_basis_);
	cdensity.fill(math::vector3<double, math::covariant>{0.0, 0.0, 0.0});

	auto iphi = 0;
	for(auto & phi : electrons.lot()){
		
		auto gphi = operations::gradient(phi);
    
    gpu::run(phi.basis().part().local_size(),
             [nst = phi.set_part().local_size(), occ = begin(electrons.occupations()[iphi]),
              ph = begin(phi.matrix()), gph = begin(gphi.matrix()), cdens = begin(cdensity.linear())] GPU_LAMBDA (auto ip){
               for(int ist = 0; ist < nst; ist++) cdens[ip] += 0.5*occ[ist]*imag(conj(ph[ip][ist])*gph[ip][ist] - conj(gph[ip][ist])*ph[ip][ist]);
             });
		iphi++;
	}
  
	return cdensity;
}

template <typename HamiltonianType>
auto current(const systems::ions & ions, systems::electrons & electrons, HamiltonianType const & ham){
  return operations::integral(current_density(ions, electrons, ham));
}

}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_OBSERVABLES_CURRENT_UNIT_TEST
#undef INQ_OBSERVABLES_CURRENT_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("observables::current", "[observables::current]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using math::vector3;

}

#endif
#endif
