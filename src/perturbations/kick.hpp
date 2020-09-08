/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__KICK
#define INQ__PERTURBATIONS__KICK

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
#include <basis/field_set.hpp>

namespace inq {
namespace perturbations {

	void kick(math::vec3d kick_field, basis::field_set<basis::real_space, complex> & phi){

		//DATAOPERATIONS LOOP 4D
		for(int ix = 0; ix < phi.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < phi.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < phi.basis().local_sizes()[2]; iz++){
					
					auto ixg = phi.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = phi.basis().cubic_dist(1).local_to_global(iy);
					auto izg = phi.basis().cubic_dist(2).local_to_global(iz);
					
					auto rr = phi.basis().rvector(ixg, iyg, izg);

					auto kick_factor = exp(complex(0.0, (kick_field|rr)));
					
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] *= kick_factor;
					}
					
				}
			}
		}
	}
}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_KICK_UNIT_TEST
#undef INQ_PERTURBATIONS_KICK_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("perturbations::kick", "[perturbations::kick]") {

	using namespace inq;
	using namespace Catch::literals;
	using math::vec3d;
	
	const int nvec = 12;

	double ecut = 31.2;

	ions::UnitCell cell(vec3d(4.2, 0.0, 0.0), vec3d(0.0, 3.5, 0.0), vec3d(0.0, 0.0, 6.4));

	basis::real_space bas(cell, input::basis::cutoff_energy(ecut));

	basis::field_set<basis::real_space, complex> phi(bas, nvec);

	perturbations::kick({0.1, 0.0, 0.0}, phi);
	
}

#endif
#endif
