/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef BASIS_FIELD
#define BASIS_FIELD

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
  class field {

  public:

		typedef type value_type;
		typedef boost::multi::array<type, basis_type::dimension> cubic_type;
		
    field(const basis_type & basis):
      cubic_(sizes(basis)),
			basis_(basis){
    }

		field(const field & coeff) = delete;
		field(field && coeff) = default;
		field & operator=(const field & coeff) = delete;
		field & operator=(field && coeff) = default;
		
		const basis_type & basis() const {
			return basis_;
		}

		const auto & cubic() const {
			return cubic_;
		}

		auto & cubic() {
			return cubic_;
		}

		const auto & operator[](const long index) const {
			return cubic_.data()[index];
		}

		auto & operator[](const long index) {
			return cubic_.data()[index];
		}

		auto data() const {
			return cubic_.data();
		}

		auto data() {
			return cubic_.data();
		}
		
	private:
		
    cubic_type cubic_;
		basis_type basis_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

TEST_CASE("Class basis::field", "[field]"){
  
}

#endif

#endif
