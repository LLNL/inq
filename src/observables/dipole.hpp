/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__DIPOLE
#define INQ__OBSERVABLES__DIPOLE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <physics/constants.hpp>

namespace inq {
namespace observables {

vector3<double> dipole(basis::field<basis::real_space, double> const & density){

	CALI_CXX_MARK_FUNCTION;
	
	vector3<double> dip = {0.0, 0.0, 0.0};

	for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
		for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
			for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
				dip += density.cubic()[ix][iy][iz]*density.basis().point_op().rvector_cartesian(ix, iy, iz);
			}
		}
	}

	if(density.basis().comm().size() > 1) density.basis().comm().all_reduce_n(dip.data(), dip.size());

	return dip*density.basis().volume_element();
	
}

vector3<double> dipole(ions::geometry const & geo, const hamiltonian::atomic_potential & atomic_pot){

	using physics::constants::proton_charge;
	
	vector3<double> dip = {0.0, 0.0, 0.0};

	for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
		auto zval = atomic_pot.pseudo_for_element(geo.atoms()[iatom]).valence_charge();
		dip += proton_charge*zval*geo.coordinates()[iatom];
	}

	return dip;
}

vector3<double> dipole(systems::ions const & ions, systems::electrons const & electrons){
	return dipole(ions.geo_, electrons.atomic_pot()) + dipole(electrons.density());
}

}
}
#endif

#ifdef INQ_OBSERVABLES_DIPOLE_UNIT_TEST
#undef INQ_OBSERVABLES_DIPOLE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b);

  basis::real_space bas(box, /*spacing =*/ 0.39770182, comm);

	basis::field<basis::real_space, double> density(bas);

	SECTION("Dipole x"){
		for(int ix = 0; ix < density.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < density.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < density.basis().local_sizes()[2]; iz++){
					
					auto rr = density.basis().point_op().rvector_cartesian(ix, iy, iz);
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
					
					auto rr = density.basis().point_op().rvector_cartesian(ix, iy, iz);
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
					
					auto rr = density.basis().point_op().rvector_cartesian(ix, iy, iz);
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
