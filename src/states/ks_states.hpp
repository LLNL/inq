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

#include <math/complex.hpp>
#include <multi/array.hpp>

namespace states {
  template <class basis_type>
  class ks_states {

  public:
    
    enum class spin_config {
      UNPOLARIZED,
      POLARIZED,
      NON_COLLINEAR
    };
        
    ks_states(const spin_config spin, const double nelectrons, const basis_type & basis){

      if(spin == spin_config::NON_COLLINEAR){
				nspinor_ = 2;
				nstates_ = ceil(nelectrons);
      } else {
				nspinor_ = 1;
				nstates_ = ceil(0.5*nelectrons);
      }

      nquantumnumbers_ = 1;
      if(spin == spin_config::POLARIZED) nquantumnumbers_ = 2;

			coeff.reextent({nquantumnumbers_, basis.rsize()[0], basis.rsize()[1], basis.rsize()[2], nstates_, nspinor_});
			
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

		boost::multi::array<complex, 6> coeff;		
		
  private:

    int nspinor_;
    int nstates_;
    int nquantumnumbers_;

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
		REQUIRE(std::get<0>(extensions(st.coeff)) == st.num_quantum_numbers());
		REQUIRE(std::get<1>(extensions(st.coeff)) == pw.rsize()[0]);
		REQUIRE(std::get<2>(extensions(st.coeff)) == pw.rsize()[1]);
		REQUIRE(std::get<3>(extensions(st.coeff)) == pw.rsize()[2]);
		REQUIRE(std::get<4>(extensions(st.coeff)) == st.num_states());
		REQUIRE(std::get<5>(extensions(st.coeff)) == st.num_spinors());
  }

  SECTION("Spin polarized"){
    
    states::ks_states<basis::plane_wave> st(states::ks_states<basis::plane_wave>::spin_config::POLARIZED, 11.0, pw);
    
    REQUIRE(st.num_spinors() == 1);
    REQUIRE(st.num_states() == 6);
    REQUIRE(st.num_quantum_numbers() == 2);
		REQUIRE(std::get<0>(extensions(st.coeff)) == st.num_quantum_numbers());
		REQUIRE(std::get<1>(extensions(st.coeff)) == pw.rsize()[0]);
		REQUIRE(std::get<2>(extensions(st.coeff)) == pw.rsize()[1]);
		REQUIRE(std::get<3>(extensions(st.coeff)) == pw.rsize()[2]);
		REQUIRE(std::get<4>(extensions(st.coeff)) == st.num_states());
		REQUIRE(std::get<5>(extensions(st.coeff)) == st.num_spinors());
  }

  SECTION("Non-collinear spin"){
    
    states::ks_states<basis::plane_wave> st(states::ks_states<basis::plane_wave>::spin_config::NON_COLLINEAR, 11.0, pw);
    
    REQUIRE(st.num_spinors() == 2);
    REQUIRE(st.num_states() == 11);
    REQUIRE(st.num_quantum_numbers() == 1);
		REQUIRE(std::get<0>(extensions(st.coeff)) == st.num_quantum_numbers());
		REQUIRE(std::get<1>(extensions(st.coeff)) == pw.rsize()[0]);
		REQUIRE(std::get<2>(extensions(st.coeff)) == pw.rsize()[1]);
		REQUIRE(std::get<3>(extensions(st.coeff)) == pw.rsize()[2]);
		REQUIRE(std::get<4>(extensions(st.coeff)) == st.num_states());
		REQUIRE(std::get<5>(extensions(st.coeff)) == st.num_spinors());
  }

  
}

#endif

#endif
