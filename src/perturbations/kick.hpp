/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__KICK
#define INQ__PERTURBATIONS__KICK

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <states/orbital_set.hpp>
#include <perturbations/gauge.hpp>
#include <perturbations/none.hpp>

namespace inq {
namespace perturbations {

class kick : public perturbations::none {

public:

	template <typename CellType>
	kick(CellType const & cell, vector3<double> const & arg_kick_field, gauge arg_gauge = gauge::mixed):
		efield_(-arg_kick_field),
		vpot_(-arg_kick_field),		
		periodicity_(cell.periodicity())
	{
		if(arg_gauge == gauge::mixed){
			for(int idir = 0; idir < periodicity_; idir++) efield_[idir] = 0.0;
			for(int idir = periodicity_; idir < 3; idir++) vpot_[idir] = 0.0;
		}

		if(arg_gauge == gauge::length){
			vpot_ = {0.0, 0.0, 0.0};
		}

		if(arg_gauge == gauge::velocity){
			efield_ = {0.0, 0.0, 0.0};
		}

		assert(efield_ + vpot_ == -arg_kick_field);
	}

	template <typename PhiType>
	void zero_step(PhiType & phi) const {

		auto cov_efield = phi.basis().cell().metric().to_covariant(efield_);
		
		gpu::run(phi.basis().local_sizes()[2], phi.basis().local_sizes()[1], phi.basis().local_sizes()[0],
						 [pop = phi.basis().point_op(), ph = begin(phi.hypercubic()), cov_efield, nst = phi.set_part().local_size()] GPU_LAMBDA (auto iz, auto iy, auto ix){
							 
							 auto rr = pop.rvector(ix, iy, iz);
							 auto kick_factor = exp(complex(0.0, dot(cov_efield, rr)));
							 for(int ist = 0; ist < nst; ist++) ph[ix][iy][iz][ist] *= kick_factor;
						 });
	}

	auto has_uniform_vector_potential() const {
		return true;
	}
	
	auto uniform_vector_potential(double /*time*/) const {
		return vpot_;
	}

private:

	vector3<double> efield_;
	vector3<double> vpot_;	
	int periodicity_;

};

}
}
#endif

#ifdef INQ_PERTURBATIONS_KICK_UNIT_TEST
#undef INQ_PERTURBATIONS_KICK_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;
	
	SECTION("finite"){
	
		const int nvec = 12;
		
		double phi_absdif = 0.0;
		double phi_dif = 0.0;
		
		parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
		
		systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).finite();
		
		CHECK(box.periodicity_value() == 0);
		
		basis::real_space bas(box, /*spacing =*/ 0.39770182, comm);
		
		CHECK(bas.cell().periodicity() == 0);
		
		basis::field_set<basis::real_space, complex> phi(bas, nvec);
		
		//Construct a field
		for(int ix = 0; ix < phi.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < phi.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < phi.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						phi.hypercubic()[ix][iy][iz][ist] = complex(cos(ist+(ix+iy+iz)), 1.3*sin(ist+(cos(ix-iy-iz))));
					}
				}
			}
		}
		
		auto phi_old = phi;
		
		auto kick = perturbations::kick(box.cell(), {0.1, 0.0, 0.0});
		
		kick.zero_step(phi);
		
		for(int ix = 0; ix < phi.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < phi.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < phi.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						phi_absdif += norm(phi.hypercubic()[ix][iy][iz][ist]) - norm(phi_old.hypercubic()[ix][iy][iz][ist]);
						phi_dif += norm(phi.hypercubic()[ix][iy][iz][ist] - phi_old.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
		CHECK(phi_absdif == Approx(0).margin(1.0e-9));
		CHECK(phi_dif > 1.0e-9);
	}
	
	SECTION("periodic"){
		systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b);
		auto kick = perturbations::kick(box.cell(), {0.1, 0.2, 0.3});

		CHECK(kick.has_uniform_vector_potential());
		CHECK(kick.uniform_vector_potential(3.0)[0] == -0.1);
		CHECK(kick.uniform_vector_potential(2.0)[1] == -0.2);
		CHECK(kick.uniform_vector_potential(1.0)[2] == -0.3);
	}

	SECTION("semi periodic"){
		systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).periodicity(2);
		auto kick = perturbations::kick(box.cell(), {0.1, 0.2, 0.3});

		CHECK(kick.has_uniform_vector_potential());
		CHECK(kick.uniform_vector_potential(3.0)[0] == -0.1);
		CHECK(kick.uniform_vector_potential(2.0)[1] == -0.2);
		CHECK(kick.uniform_vector_potential(1.0)[2] == 0.0);
	}

	SECTION("velocity gauge"){
		systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).finite();
		auto kick = perturbations::kick(box.cell(), {0.1, 0.2, 0.3}, perturbations::gauge::velocity);

		CHECK(kick.has_uniform_vector_potential());
		CHECK(kick.uniform_vector_potential(3.0)[0] == -0.1);
		CHECK(kick.uniform_vector_potential(2.0)[1] == -0.2);
		CHECK(kick.uniform_vector_potential(1.0)[2] == -0.3);
	}

	SECTION("length gauge"){
		systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).periodic();
		auto kick = perturbations::kick(box.cell(), {0.1, 0.2, 0.3}, perturbations::gauge::length);

		CHECK(kick.has_uniform_vector_potential());
		CHECK(kick.uniform_vector_potential(3.0)[0] == 0.0);
		CHECK(kick.uniform_vector_potential(2.0)[1] == 0.0);
		CHECK(kick.uniform_vector_potential(1.0)[2] == 0.0);
	}

}
#endif
