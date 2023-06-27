/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__INTERACTION
#define INPUT__INTERACTION

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <utils/merge_optional.hpp>

#include <cassert>
#include <optional>
#include <stdexcept>

#include <xc.h>

namespace inq {
namespace input {

class interaction {

public:

	enum class exchange_functional {
		NONE = 0,
		LDA = XC_LDA_X,
		PBE = XC_GGA_X_PBE,
		B = XC_GGA_X_B88,
		B3LYP = XC_HYB_GGA_XC_B3LYP,
		PBE0 = XC_HYB_GGA_XC_PBEH,
		HARTREE_FOCK = -1
	};

	enum class correlation_functional {
		NONE = 0,
		LDA_PZ = XC_LDA_C_PZ,
		PBE = XC_GGA_C_PBE,
		LYP = XC_GGA_C_LYP
	};
	
	enum class induced_vector_potential {
		NONE,
		GAUGE_FIELD
	};

private:

	std::optional<bool> hartree_potential_;
	std::optional<exchange_functional> exchange_;
	std::optional<correlation_functional> correlation_;
	std::optional<induced_vector_potential> induced_vecpot_;
	std::optional<bool> fourier_pseudo_;
	double alpha_ = 0;

public:
	
	auto non_interacting() const {
		interaction inter = *this;
		inter.hartree_potential_ = false;
		inter.exchange_ = exchange_functional::NONE;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}

	auto dft() const {
		interaction inter = *this;
		inter.hartree_potential_ = true;
		return inter;
	}
	
	auto lda() const {
		interaction inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::LDA;
		inter.correlation_ = correlation_functional::LDA_PZ;		
		return inter;
	}

	auto hartree() const {
		interaction inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::NONE;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}
	
	auto hartree_fock() const {
		interaction inter = *this;
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

	auto pbe() const {
		interaction inter = *this;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::PBE;
		inter.correlation_ = correlation_functional::PBE;
		return inter;
	}

	auto pbe0()  const {
		interaction inter = *this;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::PBE0;
		inter.correlation_ = correlation_functional::NONE;
		return inter;
	}

	auto b3lyp()  const {
		interaction inter = *this;
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

	auto real_space_pseudo() const {
		interaction inter = *this;
		inter.fourier_pseudo_ = false;
		return inter;
	}

	auto fourier_pseudo() const {
		interaction inter = *this;
		inter.fourier_pseudo_ = true;
		return inter;
	}
		
	auto fourier_pseudo_value() const {
		return fourier_pseudo_.value_or(false);
	}

	auto gauge_field() const {
		interaction inter = *this;
		inter.induced_vecpot_ = induced_vector_potential::GAUGE_FIELD;
		inter.alpha_ = -4*M_PI;
		return inter;
	}

	auto gauge_field(const double alpha){
		interaction inter = *this;
		inter.induced_vecpot_ = induced_vector_potential::GAUGE_FIELD;
		inter.alpha_ = alpha;
		return inter;
	}
	
	auto induced_vector_potential_value() const {
		return induced_vecpot_.value_or(induced_vector_potential::NONE);
	}

	auto alpha_value() const {
		return alpha_;
	}

		
};
    
}
}
#endif

#ifdef INQ_INPUT_INTERACTION_UNIT_TEST
#undef INQ_INPUT_INTERACTION_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
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

    auto inter = input::interaction{}.non_interacting().fourier_pseudo();
    
		CHECK(not inter.self_consistent());
		CHECK(inter.fourier_pseudo_value() == true);
		CHECK(inter.exchange_coefficient() == 0.0);		
  }
	
  SECTION("Hartee-Fock"){

    auto inter = input::interaction{}.hartree_fock();
		CHECK(inter.exchange_coefficient() == 1.0);
  }

}
#endif
