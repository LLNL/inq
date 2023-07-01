/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INQ__QUANTITY
#define INQ__INQ__QUANTITY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/run.hpp>
#include <math/vector3.hpp>

namespace inq {

	template <class MagnitudeType, class ElementType = double>
	class quantity {
	public:
		using magnitude = MagnitudeType;
		using element_type = ElementType;

	public:
		template<
			class Other,
			class = std::enable_if_t<not std::is_base_of<quantity, Other>::value, int>
		>
		quantity(Other const&) {static_assert(sizeof(Other*) and false, "deleted, a value is probably missing units");}

		template<class Other>
		explicit operator Other() const {static_assert(sizeof(Other*) and false, "deleted, a value is probably missing units"); return Other{};}

		quantity() = default;
		quantity(quantity const&) = default;
		quantity & operator=(quantity const&) = default;
		
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
#endif

#ifdef INQ_INQ_QUANTITY_UNIT_TEST
#undef INQ_INQ_QUANTITY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

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

