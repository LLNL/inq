/* -*- indent-tabs-mode: t -*- */

#ifndef OBSERVABLES__DIPOLE
#define OBSERVABLES__DIPOLE

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math/vec3d.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>

namespace observables {

	math::vec3d dipole(basis::field<basis::real_space, double> & density){

		math::vec3d dip = {0.0, 0.0, 0.0};
		
		//DATAOPERATIONS LOOP 3D
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto ixg = density.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = density.basis().cubic_dist(1).local_to_global(iy);
					auto izg = density.basis().cubic_dist(2).local_to_global(iz);
					
					dip += density.cubic()[ix][iy][iz]*density.basis().rvector(ixg, iyg, izg);
					
				}
			}
		}

		return dip*density.basis().volume_element();
		
	}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("observables::dipole", "[observables::dipole]") {

	using namespace Catch::literals;
	using math::vec3d;
	
  const int nvec = 12;

	double ecut = 31.2;

  ions::UnitCell cell(vec3d(4.2, 0.0, 0.0), vec3d(0.0, 3.5, 0.0), vec3d(0.0, 0.0, 6.4));

  basis::real_space bas(cell, input::basis::cutoff_energy(ecut));

	basis::field<basis::real_space, double> density(bas);

	SECTION("Dipole x"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto ixg = density.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = density.basis().cubic_dist(1).local_to_global(iy);
					auto izg = density.basis().cubic_dist(2).local_to_global(iz);
					
					auto rr = density.basis().rvector(ixg, iyg, izg);

					density.cubic()[ix][iy][iz] = rr[0]*exp(-norm(rr));						
						
				}
			}
		}
	
		auto dipole = observables::dipole(density);

		CHECK(dipole[0] == 2.6692428234_a);
		CHECK(fabs(dipole[1]) < 1e-16);
		CHECK(fabs(dipole[2]) < 1e-16);
		
	}
	
	SECTION("Dipole yz"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto ixg = density.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = density.basis().cubic_dist(1).local_to_global(iy);
					auto izg = density.basis().cubic_dist(2).local_to_global(iz);
					
					auto rr = density.basis().rvector(ixg, iyg, izg);

					density.cubic()[ix][iy][iz] = (rr[1] + 2.0*rr[2])*exp(-norm(rr));
						
				}
			}
		}
	
		auto dipole = observables::dipole(density);

		CHECK(fabs(dipole[0]) < 1e-16);
		CHECK(dipole[1] == 2.496792_a);
		CHECK(dipole[2] == 5.484786_a);
		
	}

	SECTION("Dipole xz"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto ixg = density.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = density.basis().cubic_dist(1).local_to_global(iy);
					auto izg = density.basis().cubic_dist(2).local_to_global(iz);
					
					auto rr = density.basis().rvector(ixg, iyg, izg);

					density.cubic()[ix][iy][iz] = (sin(rr[1]) + rr[0]*rr[2])*exp(-norm(rr));
						
				}
			}
		}
	
		auto dipole = observables::dipole(density);

		CHECK(dipole[0] == -0.0000688417_a);
		CHECK(dipole[1] == 2.042557045_a);
		CHECK(fabs(dipole[2]) < 1e-16);
		
	}
}

#endif
#endif
