/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__SUM
#define INQ__PERTURBATIONS__SUM

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>

namespace inq {
namespace perturbations {


template <typename PertTypeA, typename PertTypeB>
class sum {

  PertTypeA perta_;
  PertTypeB pertb_;
  
public:

  sum(PertTypeA arg_perta, PertTypeB arg_pertb):
    perta_(std::move(arg_perta)),
    pertb_(std::move(arg_pertb))
  {
  }
  
	template <typename PhiType>  
	void zero_step(PhiType & phi) const {
    perta_.zero_step(phi);
    pertb_.zero_step(phi);   
	}
	
	auto has_uniform_electric_field() const {
		return perta_.has_uniform_electric_field() or pertb_.has_uniform_electric_field();
	}

	auto uniform_electric_field(double time) const {
    vector3<double> efield({0.0, 0.0, 0.0});
    if(perta_.has_uniform_electric_field()) efield += perta_.uniform_electric_field(time);
    if(pertb_.has_uniform_electric_field()) efield += pertb_.uniform_electric_field(time);
    return efield;
	}
	
	auto has_uniform_vector_potential() const {
		return perta_.has_uniform_vector_potential() or pertb_.has_uniform_vector_potential();
	}

	auto uniform_vector_potential(double time) const {
    vector3<double> vpotential({0.0, 0.0, 0.0});
    if(perta_.has_uniform_vector_potential()) vpotential += perta_.uniform_vector_potential(time);
    if(pertb_.has_uniform_vector_potential()) vpotential += pertb_.uniform_vector_potential(time);
    return vpotential;
	}
    template<typename PotentialType>
    void potential(const double time, PotentialType & potential) const {
		perta_.potential(time, potential);
		pertb_.potential(time, potential);
    }

};

template <typename PertTypeA, typename PertTypeB>
auto operator+(PertTypeA perta, PertTypeB pertb){
  return perturbations::sum(perta, pertb);
}

}
}
#endif

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_SUM_UNIT_TEST
#undef INQ_PERTURBATIONS_SUM_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

#include <perturbations/kick.hpp>
#include <perturbations/laser.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).periodic();
  auto kick = perturbations::kick(box.cell(), {0.1, 0.2, 0.3}, perturbations::gauge::velocity);
    
  auto ps = perturbations::sum(kick, perturbations::laser({1.0, 1.0, 1.0}, 1.0_Ha));

  CHECK(ps.has_uniform_electric_field());
  CHECK(ps.uniform_electric_field(M_PI/2.0)[0] == 1.0);
  CHECK(ps.uniform_electric_field(M_PI/2.0)[1] == 1.0);
  CHECK(ps.uniform_electric_field(M_PI/2.0)[2] == 1.0);
  
  CHECK(ps.has_uniform_vector_potential());
  CHECK(ps.uniform_vector_potential(3.0)[0] == -0.1);
  CHECK(ps.uniform_vector_potential(2.0)[1] == -0.2);
  CHECK(ps.uniform_vector_potential(1.0)[2] == -0.3);

  auto ps2 = kick + kick;
  
  CHECK(not ps2.has_uniform_electric_field());

  CHECK(ps2.has_uniform_vector_potential());
  CHECK(ps2.uniform_vector_potential(3.0)[0] == -0.2);
  CHECK(ps2.uniform_vector_potential(2.0)[1] == -0.4);
  CHECK(ps2.uniform_vector_potential(1.0)[2] == -0.6);

  auto ps3 = kick + kick + perturbations::laser({1.0, 1.0, 1.0}, 1.0_Ha);

  CHECK(ps3.has_uniform_electric_field());
  CHECK(ps3.uniform_electric_field(M_PI/2.0)[0] == 1.0);
  CHECK(ps3.uniform_electric_field(M_PI/2.0)[1] == 1.0);
  CHECK(ps3.uniform_electric_field(M_PI/2.0)[2] == 1.0);
  
  CHECK(ps3.has_uniform_vector_potential());
  CHECK(ps3.uniform_vector_potential(3.0)[0] == -0.2);
  CHECK(ps3.uniform_vector_potential(2.0)[1] == -0.4);
  CHECK(ps3.uniform_vector_potential(1.0)[2] == -0.6);
  
}
#endif
