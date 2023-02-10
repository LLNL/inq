/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__INTERACTION
#define INPUT__INTERACTION

/*
	Copyright (C) 2019 Xavier Andrade

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

#include <cassert>
#include <optional>
#include <stdexcept>

#include <utils/merge_optional.hpp>

namespace inq {
namespace input {

class interaction {

public:

	// these numbers match the libxc definition
	enum class exchange_functional {
		NONE = 0,
		LDA = 1,
		PBE = 101,
		B = 106,
		B3LYP = 402,
		PBE0 = 406,
		HARTREE_FOCK = -1
	};

	enum class correlation_functional {
		NONE = 0,
		LDA_PZ = 9,
		PBE = 130,
		LYP = 131
	};


	interaction(){
	}

	static auto non_interacting(){
		interaction inter;
		inter.hartree_potential_ = false;
		inter.exchange_ = exchange_functional::NONE;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}

	static auto dft(){
		interaction inter;
		inter.hartree_potential_ = true;
		return inter;
	}
	
	static auto lda(){
		interaction inter;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::LDA;
		inter.correlation_ = correlation_functional::LDA_PZ;		
		return inter;
	}

	static auto hartree_fock(){
		interaction inter;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::HARTREE_FOCK;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}
		
	auto exchange() const {
		return exchange_.value_or(exchange_functional::PBE);
	}

	auto correlation() const {
		return correlation_.value_or(correlation_functional::PBE);
	}

	static auto pbe() {
		interaction inter;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::PBE;
		inter.correlation_ = correlation_functional::PBE;
		return inter;
	}

	static auto pbe0() {
		interaction inter;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::PBE0;
		inter.correlation_ = correlation_functional::NONE;
		return inter;
	}

	static auto b3lyp() {
		interaction inter;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::B3LYP;
		inter.correlation_ = correlation_functional::NONE;
		return inter;
	}

	auto exchange_coefficient() const {
		if(exchange() == exchange_functional::HARTREE_FOCK) return 1.0;
		if(exchange() == exchange_functional::NONE) return 0.0;
		throw std::runtime_error("inq internal error: exchange coefficient not known here for true functionals");
	}

	auto hartree_potential() const {
		return hartree_potential_.value_or(true);
	}
	
	auto self_consistent() const {
		return hartree_potential() or exchange() != exchange_functional::NONE or correlation() != correlation_functional::NONE;
	}

	static auto real_space_pseudo(){
		interaction inter;
		inter.fourier_pseudo_ = false;
		return inter;
	}

	static auto fourier_pseudo(){
		interaction inter;
		inter.fourier_pseudo_ = true;
		return inter;
	}
		
	auto fourier_pseudo_value() const {
		return fourier_pseudo_.value_or(false);
	}

	friend auto operator|(const interaction & inter1, const interaction & inter2){
		using inq::utils::merge_optional;
		
		interaction rinter;
		rinter.hartree_potential_	= merge_optional(inter1.hartree_potential_, inter2.hartree_potential_);
		rinter.exchange_	= merge_optional(inter1.exchange_, inter2.exchange_);
		rinter.correlation_	= merge_optional(inter1.correlation_, inter2.correlation_);
		rinter.fourier_pseudo_	= merge_optional(inter1.fourier_pseudo_, inter2.fourier_pseudo_);
		return rinter;
	}

private:

	std::optional<bool> hartree_potential_;
	std::optional<exchange_functional> exchange_;
	std::optional<correlation_functional> correlation_;
	std::optional<bool> fourier_pseudo_;
		
};
    
}
}
#endif

#ifdef INQ_INPUT_INTERACTION_UNIT_TEST
#undef INQ_INPUT_INTERACTION_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("class input::interaction", "[input::interaction]") {
  
	using namespace inq;
	using namespace Catch::literals;

	SECTION("Defaults"){

    input::interaction inter;

		CHECK(inter.hartree_potential() == true);
		CHECK(inter.exchange() == input::interaction::exchange_functional::PBE);
		CHECK(inter.correlation() == input::interaction::correlation_functional::PBE);
		CHECK(inter.fourier_pseudo_value() == false);
		CHECK_THROWS(inter.exchange_coefficient());
  }

  SECTION("Composition"){

    auto inter = input::interaction::non_interacting() | input::interaction::fourier_pseudo();
    
		CHECK(not inter.self_consistent());
		CHECK(inter.fourier_pseudo_value() == true);
		CHECK(inter.exchange_coefficient() == 0.0);		
  }
	
  SECTION("Hartee-Fock"){

    auto inter = input::interaction::hartree_fock();
		CHECK(inter.exchange_coefficient() == 1.0);
  }

}
#endif
