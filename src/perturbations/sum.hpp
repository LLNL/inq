/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__SUM
#define INQ__PERTURBATIONS__SUM

/*
 Copyright (C) 2023 Xavier Andrade

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
    math::vector3<double> efield({0.0, 0.0, 0.0});
    if(perta_.has_uniform_electric_field()) efield += perta_.uniform_electric_field(time);
    if(pertb_.has_uniform_electric_field()) efield += pertb_.uniform_electric_field(time);
    return efield;
	}
	
	auto has_uniform_vector_potential() const {
		return perta_.has_uniform_vector_potential() or pertb_.has_uniform_vector_potential();
	}

	auto uniform_vector_potential(double time) const {
    math::vector3<double> vpotential({0.0, 0.0, 0.0});
    if(perta_.has_uniform_vector_potential()) vpotential += perta_.uniform_vector_potential(time);
    if(pertb_.has_uniform_vector_potential()) vpotential += pertb_.uniform_vector_potential(time);
    return vpotential;
	}

};

template <typename PertTypeA, typename PertTypeB>
auto operator+(PertTypeA perta, PertTypeB pertb){
  return perturbations::sum(perta, pertb);
}

}
}

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

TEST_CASE("perturbations::sum", "[perturbations::sum]") {

  systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).periodic().cutoff_energy(10.3_Ry);
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
  
}

#endif
#endif
