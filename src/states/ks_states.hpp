#ifndef KS_STATES
#define KS_STATES

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

namespace states {
  template <class basis_type>
  class ks_states {

  public:
    
    enum class spin_config {
      UNPOLARIZED,
      POLARIZED,
      NON_COLINEAR
    };
        
    ks_states(const spin_config spin, const double nelectrons, const basis_type & basis):
      basis_(basis){

      if(spin == spin_config::NON_COLINEAR){
	nspinor_ = 2;
	nstates_ = ceil(nelectrons);
      } else {
	nspinor_ = 1;
	nstates_ = ceil(0.5*nelectrons);
      }

      nquantumnumbers_ = 1;
      if(spin == spin_config::POLARIZED) nquantumnumbers_ = 2;
      
    }

    int num_states() const {
      return nstates_;
    }
    
    int num_spinors() const {
      return nspinor_;
    }

    int num_quantum_numbers() const {
      return nquantumnumbers_;
    }
    
  private:

    int nspinor_;
    int nstates_;
    int nquantumnumbers_;
    basis_type basis_;    
    
  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/plane_wave.hpp>

TEST_CASE("Class states::ks_states", "[ks_states]"){

  using math::d3vector;
  
  double ecut = 30.0;
  double ll = 10.0;
  
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::plane_wave pw(cell, ecut);
  
  SECTION("Spin unpolarized"){
    
    states::ks_states<basis::plane_wave> st(states::ks_states<basis::plane_wave>::spin_config::UNPOLARIZED, 11.0, pw);
    
    REQUIRE(st.num_spinors() == 1);
    REQUIRE(st.num_states() == 6);
    REQUIRE(st.num_quantum_numbers() == 1);
  }

  SECTION("Spin polarized"){
    
    states::ks_states<basis::plane_wave> st(states::ks_states<basis::plane_wave>::spin_config::POLARIZED, 11.0, pw);
    
    REQUIRE(st.num_spinors() == 1);
    REQUIRE(st.num_states() == 6);
    REQUIRE(st.num_quantum_numbers() == 2);
  }

  SECTION("Non-colinear spin"){
    
    states::ks_states<basis::plane_wave> st(states::ks_states<basis::plane_wave>::spin_config::NON_COLINEAR, 11.0, pw);
    
    REQUIRE(st.num_spinors() == 2);
    REQUIRE(st.num_states() == 11);
    REQUIRE(st.num_quantum_numbers() == 1);
  }

  
}

#endif

#endif
