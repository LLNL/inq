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

	template <class type1, unsigned long dim, class type2>
	static auto array_append(const std::array<type1, dim> & array, type2 element){
    std::array<type1, dim + 1> result;
    std::copy(array.cbegin(), array.cend(), result.begin());
		result[dim] = element;
    return result;
	}
	
	template<class basis_type, class type>
  class coefficients_set {

  public:

		typedef type value_type;
		typedef boost::multi::array<type, basis_type::dimension + 1> cubic_type;
		
    coefficients_set(const basis_type & basis, const int num_vectors):
			cubic_(array_append(sizes(basis), num_vectors)),
			num_vectors_(num_vectors),
			basis_(basis){
    }

		coefficients_set(const coefficients_set & coeff) = delete;
		coefficients_set(coefficients_set && coeff) = default;
		coefficients_set & operator=(const coefficients_set & coeff) = delete;
		coefficients_set & operator=(coefficients_set && coeff) = default;
		
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

TEST_CASE("Class basis::coefficients_set", "[coefficients_set]"){
  
}

#endif

#endif
