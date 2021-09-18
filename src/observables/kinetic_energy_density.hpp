/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__KINETIC_ENERGY_DENSITY
#define INQ__OBSERVABLES__KINETIC_ENERGY_DENSITY

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

#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <systems/electrons.hpp>
#include <operations/gradient.hpp>

namespace inq {
namespace observables {

basis::field<basis::real_space, double> kinetic_energy_density(systems::electrons const & electrons){

	CALI_CXX_MARK_FUNCTION;

  auto gphi = operations::gradient(electrons.phi_.fields());
	basis::field<basis::real_space, double> density(electrons.phi_.fields().basis());

	gpu::run(density.basis().part().local_size(),
					 [nst = gphi.set_part().local_size(),
						occ = begin(electrons.phi_.occupations()),
						gph = begin(gphi.matrix()),
						den = begin(density.linear())]
					 GPU_LAMBDA (auto ipoint){
						 den[ipoint] = 0.0;
						 for(int ist = 0; ist < nst; ist++) den[ipoint] += occ[ist]*norm(gph[ipoint][ist]);
						 den[ipoint] *= 0.5;
					 });
	
	if(gphi.set_comm().size() > 1){
		CALI_CXX_MARK_SCOPE("kinetic_energy_density::reduce");
		gphi.set_comm().all_reduce_in_place_n(raw_pointer_cast(density.linear().data_elements()), density.linear().size(), std::plus<>{});
	}

	return density;
	
}

}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_OBSERVABLES_KINETIC_ENERGY_DENSITY_UNIT_TEST
#undef INQ_OBSERVABLES_KINETIC_ENERGY_DENSITY_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("observables::kinetic_energy_density", "[observables::kinetic_energy_density]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;

}

#endif
#endif
