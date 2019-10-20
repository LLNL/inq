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
                                        PBE = 30,
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

		auto exchange_coefficient() const {
			if(theory_ == electronic_theory::HARTREE_FOCK) return 1.0;
			return 0.0;
		}
		
  private:

    nonstd::optional<electronic_theory> theory_;
    nonstd::optional<exchange_functional> exchange_;
    nonstd::optional<correlation_functional> correlation_;
    
  };
    
}

////////////////////////////////////////////////////////

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#endif
   
#endif
