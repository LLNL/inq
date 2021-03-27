/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INQ__QUANTITY
#define INQ__INQ__QUANTITY

/*
 Copyright (C) 2021 Xavier Andrade

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

namespace inq {

	template <class MagnitudeType, class ElementType = double>

	class quantity {
	public:
		using magnitude = MagnitudeType;
		using element_type = ElementType;

		static auto from_atomic_units(element_type const & au_value){
			quantity qq;
			qq.value_ = au_value;
			return qq;
		}

		static auto zero(){
			return from_atomic_units(0.0);
		}

		auto in_atomic_units() const {
			return value_;
		}

		friend auto operator*(double scal, quantity quant){
			quant.value_ *= scal;
			return quant;
		}

		friend auto operator*(quantity quant, double scal){
			quant.value_ *= scal;
			return quant;
		}

	private:

		element_type value_;
		
	};
	
}

#ifdef INQ_INQ_QUANTITY_UNIT_TEST
#undef INQ_INQ_QUANTITY_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("inq::quantity", "[inq::quantity]") {

	using namespace inq;
	using namespace Catch::literals;

	auto zz = quantity<void>::zero();
	CHECK(zz.in_atomic_units() == 0.0_a);

	auto rr = quantity<void>::from_atomic_units(25.5);
	CHECK(rr.in_atomic_units() == 25.5_a);

	auto rr2 = 4.0*rr;
	CHECK(rr2.in_atomic_units() == 102.0_a);

	rr2 = rr*4.0;
	CHECK(rr2.in_atomic_units() == 102.0_a);
}

#endif

#endif

