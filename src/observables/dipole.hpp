/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__DIPOLE
#define INQ__OBSERVABLES__DIPOLE

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
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <physics/constants.hpp>

namespace inq {
namespace observables {

math::vector3<double> dipole(basis::field<basis::real_space, double> & density){
	
	math::vector3<double> dip = {0.0, 0.0, 0.0};
	
	//DATAOPERATIONS LOOP 3D
	for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
		for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
			for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
				dip += density.cubic()[ix][iy][iz]*density.basis().rvector(ix, iy, iz);
			}
		}
	}
	
	density.basis().comm().all_reduce_in_place_n(dip.data(), dip.size(), std::plus<>{});
	
	return dip*density.basis().volume_element();
	
}

math::vector3<double> dipole(ions::geometry const & geo, const hamiltonian::atomic_potential & atomic_pot){

	using physics::constants::proton_charge;
	
	math::vector3<double> dip = {0.0, 0.0, 0.0};

	for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
		auto zval = atomic_pot.pseudo_for_element(geo.atoms()[iatom]).valence_charge();
		dip += proton_charge*zval*geo.coordinates()[iatom];
	}

	return dip;
}

math::vector3<double> dipole(systems::ions const & ions, systems::electrons & electrons){
	return dipole(ions.geo_, electrons.atomic_pot_) + dipole(electrons.density_);
}


}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_OBSERVABLES_DIPOLE_UNIT_TEST
#undef INQ_OBSERVABLES_DIPOLE_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("observables::dipole", "[observables::dipole]") {

	using namespace inq;
	using namespace Catch::literals;
	using math::vector3;
	
	double ecut = 31.2;

  ions::UnitCell cell(vector3<double>(4.2, 0.0, 0.0), vector3<double>(0.0, 3.5, 0.0), vector3<double>(0.0, 0.0, 6.4));

  basis::real_space bas(cell, input::basis::cutoff_energy(ecut));

	basis::field<basis::real_space, double> density(bas);

	SECTION("Dipole x"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto rr = density.basis().rvector(ix, iy, iz);
					density.cubic()[ix][iy][iz] = rr[0]*exp(-norm(rr));
						
				}
			}
		}
	
		auto dipole = observables::dipole(density);

		CHECK(dipole[0] == 2.6692428234_a);
		CHECK(fabs(dipole[1]) < 1e-14);
		CHECK(fabs(dipole[2]) < 1e-14);
		
	}
	
	SECTION("Dipole yz"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto rr = density.basis().rvector(ix, iy, iz);
					density.cubic()[ix][iy][iz] = (rr[1] + 2.0*rr[2])*exp(-norm(rr));
						
				}
			}
		}
	
		auto dipole = observables::dipole(density);

		CHECK(fabs(dipole[0]) < 1e-14);
		CHECK(dipole[1] == 2.496792_a);
		CHECK(dipole[2] == 5.484786_a);
		
	}

	SECTION("Dipole xz"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto rr = density.basis().rvector(ix, iy, iz);
					density.cubic()[ix][iy][iz] = (sin(rr[1]) + rr[0]*rr[2])*exp(-norm(rr));
						
				}
			}
		}
	
		auto dipole = observables::dipole(density);

		CHECK(dipole[0] == -0.0000688417_a);
		CHECK(dipole[1] == 2.042557045_a);
		CHECK(fabs(dipole[2]) < 1e-14);
		
	}
}

#endif
#endif
