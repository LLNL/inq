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

#include <gpu/run.hpp>
#include <math/vector3.hpp>

namespace inq {

	template <class MagnitudeType, class ElementType = double>
	class quantity {
	public:
		using magnitude = MagnitudeType;
		using element_type = ElementType;

	public:
		
		GPU_FUNCTION static auto from_atomic_units(element_type const & au_value){
			quantity qq;
			qq.value_ = au_value;
			return qq;
		}

		GPU_FUNCTION static auto zero(){
			return from_atomic_units(0.0);
		}

		GPU_FUNCTION auto in_atomic_units() const {
			return value_;
		}
		
		GPU_FUNCTION auto operator*=(double scal){
			value_ *= scal;
			return *this;
		}
		
		GPU_FUNCTION friend auto operator*(double scal, quantity quant){
			quant.value_ *= scal;
			return quant;
		}

		GPU_FUNCTION friend auto operator*(quantity quant, double scal){
			quant.value_ *= scal;
			return quant;
		}

		GPU_FUNCTION friend auto operator/(quantity quant, double scal){
			quant.value_ /= scal;
			return quant;
		}

		GPU_FUNCTION friend auto operator/(quantity num, quantity den){
			return num.value_/den.value_;
		}

		GPU_FUNCTION auto operator-() const {
			return from_atomic_units(-value_);
		}

		GPU_FUNCTION auto operator+=(quantity quant){
			value_ += quant.value_;
			return *this;
		}
				
		GPU_FUNCTION auto operator+(quantity quant) const {
			quant += *this;
			return quant;
		}
		
		GPU_FUNCTION auto operator-=(quantity quant){
			value_ -= quant.value_;
			return *this;
		}
				
		GPU_FUNCTION auto operator-(quantity quant) const {
			quant -= *this;
			return quant;
		}

		GPU_FUNCTION auto operator==(quantity quant) const {
			return value_ == quant.value_;
		}
		
		GPU_FUNCTION auto operator<=(quantity quant) const {
			return value_ <= quant.value_;
		}

		GPU_FUNCTION auto operator>=(quantity quant) const {
			return value_ >= quant.value_;
		}

		GPU_FUNCTION auto operator<(quantity quant) const {
			return value_ < quant.value_;
		}

		GPU_FUNCTION auto operator>(quantity quant) const {
			return value_ > quant.value_;
		}
		
		friend std::ostream& operator <<(std::ostream & out, quantity const & qq){
			out << qq.in_atomic_units();
			return out;
		}
		
	private:

		element_type value_;
		
	};

template <class MagnitudeType>
GPU_FUNCTION auto in_atomic_units(quantity<MagnitudeType> const & qua) {
	return qua.in_atomic_units();
}

GPU_FUNCTION auto in_atomic_units(double const & value) {
	return value;
}

GPU_FUNCTION auto in_atomic_units(vector3<double> const & value) {
	return value;
}

template <typename MagnitudeType>
GPU_FUNCTION auto in_atomic_units(vector3<MagnitudeType> const & quant) {
	return vector3<double>{quant[0].in_atomic_units(), quant[1].in_atomic_units(), quant[2].in_atomic_units()};
}

}

#ifdef INQ_INQ_QUANTITY_UNIT_TEST
#undef INQ_INQ_QUANTITY_UNIT_TEST

#include <catch2/catch_all.hpp>

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

	rr += rr2;
	CHECK(rr.in_atomic_units() == 127.5_a);

	auto rr3 = rr + rr2;
	CHECK(rr3.in_atomic_units() == 229.5_a);

	rr3 -= rr2;
	CHECK(rr3 == rr);

	rr3 = -rr2;
	CHECK(rr3.in_atomic_units() == -102.0_a);	

	rr3 *= -0.4;
	CHECK(rr3.in_atomic_units() == 40.8_a);	

	auto rr4 = rr - rr;
	CHECK(rr4.in_atomic_units() == 0.0_a);

	CHECK(in_atomic_units(rr4) == 0.0_a);

	CHECK(in_atomic_units(10.0) == 10.0_a);

	CHECK(in_atomic_units({10.0, -2.1, 4.6}) == vector3<double>{10.0, -2.1, 4.6});

	CHECK(in_atomic_units(vector3<decltype(rr)>{rr, rr, rr}) == vector3<double>{127.5, 127.5, 127.5});

}

#endif

#endif

