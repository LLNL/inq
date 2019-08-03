/* -*- indent-tabs-mode: t -*- */

#ifndef BASIS_FIELD_SET
#define BASIS_FIELD_SET

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

	template <class type1, unsigned long dim, class type2>
	static auto array_append(const std::array<type1, dim> & array, type2 element){
    std::array<type1, dim + 1> result;
    std::copy(array.cbegin(), array.cend(), result.begin());
		result[dim] = element;
    return result;
	}
	
	template<class basis_type, class type>
  class field_set {

  public:

		typedef type value_type;
		typedef boost::multi::array<type, basis_type::dimension + 1> cubic_type;
		
    field_set(const basis_type & basis, const int num_vectors):
			cubic_(array_append(sizes(basis), num_vectors)),
			num_vectors_(num_vectors),
			basis_(basis){
    }

		field_set(const field_set & coeff) = delete;
		field_set(field_set && coeff) = default;
		field_set & operator=(const field_set & coeff) = delete;
		field_set & operator=(field_set && coeff) = default;
		
		const basis_type & basis() const {
			return basis_;
		}

		const int & set_size() const {
			return num_vectors_;
		}
		
		const auto & cubic() const {
			return cubic_;
		}

		auto & cubic() {
			return cubic_;
		}

		auto operator[](long index) const {
			return &cubic_.data()[num_vectors_*index];
		}

		auto operator[](long index) {
			return &cubic_.data()[num_vectors_*index];
		}	

		auto data() const {
			return cubic_.data();
		}

		auto data() {
			return cubic_.data();
		}
		
	private:

    cubic_type cubic_;
		int num_vectors_;
		basis_type basis_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

TEST_CASE("Class basis::field_set", "[field_set]"){
  
}

#endif

#endif
