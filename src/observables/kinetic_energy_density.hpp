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

	basis::field<basis::real_space, double> density(electrons.states_basis());

	density.fill(0.0);
	
	for(auto & phi : electrons.lot()){
		auto gphi = operations::gradient(phi);
		
		gpu::run(density.basis().part().local_size(),
						 [nst = gphi.set_part().local_size(),
							occ = begin(electrons.occupations()[0]),
							gph = begin(gphi.matrix()),
							den = begin(density.linear()),
							metric = density.basis().cell().metric()]
						 GPU_LAMBDA (auto ipoint){
							 for(int ist = 0; ist < nst; ist++) den[ipoint] += 0.5*occ[ist]*metric.norm(gph[ipoint][ist]);
						 });

	}

	density.all_reduce(electrons.lot_states_comm());
	return density;
	
}

}
}
#endif

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_OBSERVABLES_KINETIC_ENERGY_DENSITY_UNIT_TEST
#undef INQ_OBSERVABLES_KINETIC_ENERGY_DENSITY_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;

}
#endif
