/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef BASIS_COEFFICIENTS_SET
#define BASIS_COEFFICIENTS_SET

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

#include <multi/array.hpp>

namespace basis {

	template<class basis_type, class type>
  class coefficients_set {

  public:

    coefficients_set(const basis_type & basis, const int num_vectors):
			num_vectors_(num_vectors),
			basis_(basis),
      cubic({basis.rsize()[0], basis.rsize()[1], basis.rsize()[2], num_vectors_}),
      linear(cubic.data(), {basis.size(), num_vectors_}){
    }

		int num_vectors_;
		const basis_type & basis_;
    boost::multi::array<type, 4> cubic;
    boost::multi::array_ref<type, 2>  linear;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

TEST_CASE("Class basis::coefficients_set", "[coefficients_set]"){
  
}

#endif

#endif
