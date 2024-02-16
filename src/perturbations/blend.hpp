/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__BLEND
#define INQ__PERTURBATIONS__BLEND

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <perturbations/absorbing.hpp>
#include <perturbations/kick.hpp>
#include <perturbations/laser.hpp>
#include <perturbations/none.hpp>

namespace inq {
namespace perturbations {

class blend {

  using any = std::variant<absorbing, kick, laser, none>;
  std::vector<any> perts_;
  
public:

  blend() = default;

  void clear() {
    perts_.clear();
  }

  template <class PertType>
  void add(PertType && pert){
    perts_.emplace_back(std::forward<PertType>(pert));
  }

  template <class PertType>
  void add(PertType const & pert){
    perts_.emplace_back(pert);
  }
	
	auto size() const {
		return (long) perts_.size();
	}

	auto has_uniform_electric_field() const {
		for(auto & pert : perts_){
			auto has = std::visit([&](auto per) { return per.has_uniform_electric_field(); }, pert);
			if(has) return true;
		}
		return false;
	}

	auto uniform_electric_field(double time) const {
    auto total = vector3<double>{0.0, 0.0, 0.0};
		for(auto & pert : perts_){
			if(std::visit([&](auto per) { return per.has_uniform_electric_field(); }, pert)){
				total += std::visit([&](auto per) { return per.uniform_electric_field(time);}, pert);
			}
		}
    return total;
	}

	auto has_uniform_vector_potential() const {
		for(auto & pert : perts_){
			auto has = std::visit([&](auto per) { return per.has_uniform_vector_potential(); }, pert);
			if(has) return true;
		}
		return false;
	}

	auto uniform_vector_potential(double time) const {
    auto total = vector3<double>{0.0, 0.0, 0.0};
		for(auto & pert : perts_){
			if(std::visit([&](auto per) { return per.has_uniform_vector_potential(); }, pert)){
				total += std::visit([&](auto per) { return per.uniform_vector_potential(time);}, pert);
			}
		}
    return total;
	}
};
	
}
}
#endif

#ifdef INQ_PERTURBATIONS_BLEND_UNIT_TEST
#undef INQ_PERTURBATIONS_BLEND_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	perturbations::blend manyp;

	auto cell = systems::cell::orthorhombic(4.2_b, 3.5_b, 6.4_b).periodic();
  auto kick = perturbations::kick(cell, {0.1, 0.2, 0.3}, perturbations::gauge::velocity);

	SECTION("2 elements"){
		auto ps = perturbations::blend{};
		ps.add(kick);
		ps.add(perturbations::laser({1.0, 1.0, 1.0}, 1.0_Ha));

		CHECK(ps.size() == 2);
		CHECK(ps.has_uniform_electric_field());

		CHECK(ps.has_uniform_electric_field());
		CHECK(ps.uniform_electric_field(M_PI/2.0)[0] == 1.0);
		CHECK(ps.uniform_electric_field(M_PI/2.0)[1] == 1.0);
		CHECK(ps.uniform_electric_field(M_PI/2.0)[2] == 1.0);

		CHECK(ps.has_uniform_vector_potential());
		CHECK(ps.uniform_vector_potential(3.0)[0] == -0.1);
		CHECK(ps.uniform_vector_potential(2.0)[1] == -0.2);
		CHECK(ps.uniform_vector_potential(1.0)[2] == -0.3);
		
	}

	SECTION("repeated elements"){
		auto ps = perturbations::blend{};
		ps.add(kick);
		ps.add(kick);
		
		CHECK(ps.size() == 2);
		CHECK(not ps.has_uniform_electric_field());

		CHECK(ps.has_uniform_vector_potential());
		CHECK(ps.uniform_vector_potential(3.0)[0] == -0.2);
		CHECK(ps.uniform_vector_potential(2.0)[1] == -0.4);
		CHECK(ps.uniform_vector_potential(1.0)[2] == -0.6);
	}

	SECTION("3 elements"){
		auto ps = perturbations::blend{};
		ps.add(kick);
		ps.add(kick);
		ps.add(perturbations::laser({1.0, 1.0, 1.0}, 1.0_Ha));

		CHECK(ps.size() == 3);

		CHECK(ps.has_uniform_electric_field());
		CHECK(ps.uniform_electric_field(M_PI/2.0)[0] == 1.0);
		CHECK(ps.uniform_electric_field(M_PI/2.0)[1] == 1.0);
		CHECK(ps.uniform_electric_field(M_PI/2.0)[2] == 1.0);
		  
		CHECK(ps.has_uniform_vector_potential());
		CHECK(ps.uniform_vector_potential(3.0)[0] == -0.2);
		CHECK(ps.uniform_vector_potential(2.0)[1] == -0.4);
		CHECK(ps.uniform_vector_potential(1.0)[2] == -0.6);
 
	}
	
}
#endif
