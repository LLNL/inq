/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__CONFIG
#define INQ__INPUT__CONFIG

/*
 Copyright (C) 2019-2023 Xavier Andrade, Alfredo A. Correa

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

#include <magnitude/energy.hpp>
#include <utils/merge_optional.hpp>

#include <cassert>
#include <optional>

namespace inq {
namespace input {

class config {

public:
	
	static auto extra_states(int value){
		config conf;
		conf.extra_states_ = value;
		return conf;
	}
	
	auto extra_states_val() const {
		return extra_states_.value_or(0);
	}

	static auto excess_charge(double value){
		config conf;
		conf.excess_charge_ = value;
		return conf;
	}

	auto excess_charge_val() const {
		return excess_charge_.value_or(0.0);
	}

	static auto temperature(quantity<magnitude::energy> value){
		config conf;
		conf.temperature_ = value;
		return conf;
	}

	auto temperature_val() const {
		return temperature_.value_or(quantity<magnitude::energy>::zero()).in_atomic_units();
	}
	
	friend auto operator|(config const & conf1, config const & conf2){
		using inq::utils::merge_optional;
		
		config rconf;
		rconf.extra_states_	= merge_optional(conf1.extra_states_, conf2.extra_states_);
		rconf.excess_charge_	= merge_optional(conf1.excess_charge_, conf2.excess_charge_);
		rconf.temperature_	= merge_optional(conf1.temperature_, conf2.temperature_);
		return rconf;
	}
	
private:
	
	std::optional<int> extra_states_;
	std::optional<double> excess_charge_;
	std::optional<quantity<magnitude::energy>> temperature_;

};

}
}

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_CONFIG_UNIT_TEST
#undef INQ_INPUT_CONFIG_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("input::config", "[input::config]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
