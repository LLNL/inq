#ifndef IONS_INTERACTION
#define IONS_INTERACTION

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

#include 

namespace states {
  template <class basis_type>
  class ks_states {

    enum class spin_config {
      UNPOLARIZED,
      POLARIZED,
      NON_COLLINEAR
    };
    
    
    ks_states(const sping_config spin, const double nelectrons, const basis_type & basis):
      basis_(basis){

      if(spin == sping_config::NON_COLLINEAR){
	nstates_ = ceil(nelectrons);
      } else {
	nstates_ = ceil(0.5*nelectrons);
      }

      nquantumnumbers_ = 1;
      if(spin == sping_config::POLARIZED) nquantumnumbers_ = 2;
      
    }
    
  private:

    int nstates_;
    int nquantumnumbers_;
    basis_type basis_;    
    
  };

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/planewave.hpp>

TEST_CASE("Function states::ks_states", "[ks_states]") {

  double ecut = 30.0;
  double ll = 10.0;
    
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::plane_wave pw(cell, ecut);

  states::ks_states<basis::plane_wave> st(states::ks_states::spin_config::UNPOLARIZED, 10.0, pw);
    
}

#endif

#endif
