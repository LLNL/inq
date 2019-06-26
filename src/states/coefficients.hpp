/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef STATES_COEFFICIENTS
#define STATES_COEFFICIENTS

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
#include <basis/real_space.hpp>
#include <states/ks_states.hpp>

namespace states {
  class coefficients {

  public:

    typedef ks_states::coeff_type type;

    coefficients(const ks_states & st, const basis::real_space & basis):
      cubic(st.cubic_dims(basis.rsize())),
      linear(cubic.data(), st.linear_dims(basis.rsize())){

    }
    
    boost::multi::array<ks_states::coeff_type, 4> cubic;
    boost::multi::array_ref<ks_states::coeff_type, 2>  linear;
    
  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

TEST_CASE("Class states::coefficients", "[coefficients]"){
  
}

#endif

#endif
