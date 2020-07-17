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

#include <nonstd/optional.hpp>
#include <cassert>

namespace inq {
namespace input {

  class interaction {

  public:

    enum class electronic_theory { NON_INTERACTING,
                                   DENSITY_FUNCTIONAL,
																	 HARTREE_FOCK
    };

    // these numbers match the libxc definition
    enum class exchange_functional { LDA = 1,
                                     PBE = 101,
                                     B = 106
    };

    enum class correlation_functional { LDA_PZ = 9,
                                        PBE = 130,
                                        LYP = 131
    };

    
    interaction(){
		}

    static auto theory(electronic_theory arg_theory){
      interaction inter;
      inter.theory_ = arg_theory;
      return inter;
    }

    static auto non_interacting(){
      interaction inter;
      inter.theory_ = electronic_theory::NON_INTERACTING;
      return inter;
    }

    static auto dft(){
      interaction inter;
      inter.theory_ = electronic_theory::DENSITY_FUNCTIONAL;
      return inter;
    }

		static auto hartree_fock(){
      interaction inter;
      inter.theory_ = electronic_theory::HARTREE_FOCK;
      return inter;
    }
		
    auto theory() const {
      return theory_.value_or(electronic_theory::DENSITY_FUNCTIONAL);
    }

    auto exchange() const {
      return exchange_.value_or(exchange_functional::LDA);
    }

    auto correlation() const {
      return correlation_.value_or(correlation_functional::LDA_PZ);
    }

		static auto pbe() {
			interaction inter;
			inter.exchange_ = exchange_functional::PBE;
			inter.correlation_ = correlation_functional::PBE;
			return inter;
		}

		auto exchange_coefficient() const {
			if(theory_ == electronic_theory::HARTREE_FOCK) return 1.0;
			return 0.0;
		}

		auto self_consistent() const {
			return theory_ != electronic_theory::NON_INTERACTING;
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
			rinter.theory_	= merge_optional(inter1.theory_, inter2.theory_);
			rinter.exchange_	= merge_optional(inter1.exchange_, inter2.exchange_);
			rinter.correlation_	= merge_optional(inter1.correlation_, inter2.correlation_);
			rinter.fourier_pseudo_	= merge_optional(inter1.fourier_pseudo_, inter2.fourier_pseudo_);
			return rinter;
		}
    
  private:

    nonstd::optional<electronic_theory> theory_;
    nonstd::optional<exchange_functional> exchange_;
    nonstd::optional<correlation_functional> correlation_;
		nonstd::optional<bool> fourier_pseudo_;
		
  };
    
}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class input::interaction", "[input::interaction]") {
  
	using namespace inq;
	using namespace Catch::literals;

	SECTION("Defaults"){

    input::interaction inter;

    CHECK(inter.theory() == input::interaction::electronic_theory::DENSITY_FUNCTIONAL);
		CHECK(inter.exchange() == input::interaction::exchange_functional::LDA);
		CHECK(inter.correlation() == input::interaction::correlation_functional::LDA_PZ);
    CHECK(inter.fourier_pseudo_value() == false);
  }

  SECTION("Composition"){

    auto inter = input::interaction::non_interacting() | input::interaction::fourier_pseudo();
    
		CHECK(not inter.self_consistent());
		CHECK(inter.fourier_pseudo_value() == true);
  }

}

#endif
   
#endif
