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

basis::field<basis::real_space, vector3<double>> dipole_density(basis::field<basis::real_space, double> const & density) {
	CALI_CXX_MARK_FUNCTION;

	basis::field<basis::real_space, vector3<double>> dipole(density.basis());

	gpu::run(density.basis().local_sizes()[2], density.basis().local_sizes()[1], density.basis().local_sizes()[0],
					 [point_op = density.basis().point_op(), dens = begin(density.cubic()), dip = begin(dipole.cubic())] GPU_LAMBDA (auto iz, auto iy, auto ix){

						 dip[ix][iy][iz] = point_op.rvector_cartesian(ix, iy, iz)*dens[ix][iy][iz];
					 });
	
	return dipole;
}

auto dipole(basis::field<basis::real_space, double> const & density) {
	CALI_CXX_MARK_FUNCTION;
	return operations::integral(dipole_density(density));
}

vector3<double> dipole(systems::ions const & ions, const hamiltonian::atomic_potential & atomic_pot){

	using physics::constants::proton_charge;
	
	vector3<double> dip = {0.0, 0.0, 0.0};

	for(int iatom = 0; iatom < ions.size(); iatom++){
		auto zval = atomic_pot.pseudo_for_element(ions.atoms()[iatom]).valence_charge();
		dip += proton_charge*zval*ions.positions()[iatom];
	}

	return dip;
}

vector3<double> dipole(systems::ions const & ions, systems::electrons const & electrons){
	return dipole(ions, electrons.atomic_pot()) + dipole(electrons.density());
}

}
}
#endif

#ifdef INQ_OBSERVABLES_DIPOLE_UNIT_TEST
#undef INQ_OBSERVABLES_DIPOLE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	if(comm.size() > 4) return;
	
  basis::real_space bas(systems::cell::orthorhombic(4.2_b, 3.5_b, 6.4_b), /*spacing =*/ 0.39770182, comm);

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
