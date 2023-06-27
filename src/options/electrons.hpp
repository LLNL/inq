/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPTIONS__ELECTRONS
#define INQ__OPTIONS__ELECTRONS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <magnitude/energy.hpp>
#include <states/ks_states.hpp>
#include <utils/merge_optional.hpp>

#include <cassert>
#include <optional>

namespace inq {
namespace options {

class electrons {

	std::optional<int> extra_states_;
	std::optional<double> excess_charge_;
	std::optional<quantity<magnitude::energy>> temperature_;
	std::optional<states::ks_states::spin_config> spin_;
	std::optional<double> spacing_;
	std::optional<bool> double_grid_;	
	std::optional<double> density_factor_;
	std::optional<bool> spherical_grid_;

public:
	
	auto extra_states(int value){
		electrons conf = *this;
		conf.extra_states_ = value;
		return conf;
	}
	
	auto extra_states_val() const {
		return extra_states_.value_or(0);
	}

	auto excess_charge(double value){
		electrons conf = *this;
		conf.excess_charge_ = value;
		return conf;
	}

	auto excess_charge_val() const {
		return excess_charge_.value_or(0.0);
	}

	auto temperature(quantity<magnitude::energy> value){
		electrons conf = *this;
		conf.temperature_ = value;
		return conf;
	}

	auto temperature_val() const {
		return temperature_.value_or(quantity<magnitude::energy>::zero()).in_atomic_units();
	}

	auto spin_unpolarized(){
		electrons conf = *this;
		conf.spin_ = states::ks_states::spin_config::UNPOLARIZED;
		return conf;
	}
	
	auto spin_polarized(){
		electrons conf = *this;
		conf.spin_ = states::ks_states::spin_config::POLARIZED;
		return conf;
	}

	auto spin_orbit(){
		electrons conf = *this;
		conf.spin_ = states::ks_states::spin_config::NON_COLLINEAR;
		return conf;
	}
	
	auto spin_val() const {
		return spin_.value_or(states::ks_states::spin_config::UNPOLARIZED);
	}
	
	auto num_spin_components_val() const {
		if(spin_val() == states::ks_states::spin_config::POLARIZED) return 2;
		return 1;
	}

	auto cutoff(quantity<magnitude::energy> arg_ecut){
		electrons conf = *this;
		conf.spacing_ = M_PI*sqrt(0.5/arg_ecut.in_atomic_units());
		return conf;		
	}

	auto spacing(quantity<magnitude::length> arg_spacing){
		electrons conf = *this;
		conf.spacing_ = arg_spacing.in_atomic_units();
		return conf;		
	}

	auto spacing_value() const {
		if(not spacing_.has_value()) throw std::runtime_error("Error: the cutoff energy or the spacing have not been set");
		return *spacing_;
	}

	auto double_grid(){
		electrons conf = *this;
		conf.double_grid_ = true;
		return conf;				
	}
	
	auto double_grid_value() const {
		return double_grid_.value_or(false);
	}

	auto density_factor(double arg_factor){
		electrons conf = *this;
		conf.density_factor_ = arg_factor;
		return conf;
	}

	auto density_factor_value() const {
		return density_factor_.value_or(1.0);
	}

	auto spherical_grid(bool arg_sph_grid){
		electrons conf = *this;		
		conf.spherical_grid_ = arg_sph_grid;
		return conf;
	}
	
	auto spherical_grid_value() const {
		return spherical_grid_.value_or(false);
	}

};

}
}
#endif

#ifdef INQ_OPTIONS_ELECTRONS_UNIT_TEST
#undef INQ_OPTIONS_ELECTRONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;

	auto conf = options::electrons{}.spacing(23.0_b);
	CHECK(conf.spacing_value() == 23.0_a);

	
}
#endif
