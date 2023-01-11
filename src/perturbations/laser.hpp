/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__LASER
#define INQ__PERTURBATIONS__LASER

/*
 Copyright (C) 2019-2023 Xavier Andrade, Alfredo Correa, Yifan Yao.

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

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>
#include <perturbations/gauge.hpp>

namespace inq {
namespace perturbations {

class laser {

	math::vector3<double, math::cartesian> polarization_;
	double frequency_;
	
public:
	laser(math::vector3<double, math::cartesian> polarization, quantity<magnitude::energy> frequency,gauge arg_gauge = gauge::length):
		polarization_(polarization),
		frequency_(frequency.in_atomic_units()),
		gauge_(arg_gauge) {
		assert(gauge_ != gauge::mixed);
	}
	
	template <typename DummyType>
	void zero_step(DummyType &) const {
	}
	
	auto has_uniform_electric_field() const {
		return gauge_ == gauge::length;
	}
	
	auto uniform_electric_field(double time) const {
		return polarization_*sin(time*frequency_);
	}
	
	auto has_uniform_vector_potential() const {
		return gauge_ == gauge::velocity;
	}
	
	auto uniform_vector_potential(double time) const {
		//E=-1/c*dA/dt
		return polarization_/frequency_*(cos(time*frequency_) - 1.0);
	}
	
	template <typename OutputStream>
	void print_info(OutputStream & out) {
		auto freq_ev = frequency_*27.211383;
		
		out << "Frequency :    " << frequency_ << " Ha" << std::endl;
		out << "               " << freq_ev << " eV" << std::endl;
		out << "               " << freq_ev*241.7991 << " THz" << std::endl;
		out << "               " << 1239.84193/freq_ev << " nm" << std::endl;
		
	}
private:
	gauge gauge_;

};

}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_LASER_UNIT_TEST
#undef INQ_PERTURBATIONS_LASER_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

using namespace inq;
using namespace Catch::literals;
using namespace magnitude;
using math::vector3;

TEST_CASE("perturbations::laser", "[perturbations::laser]") {

	perturbations::laser las({1.0, 0.0, 0.0}, 1.0_eV);

	las.print_info(std::cout);
	
	CHECK(las.has_uniform_electric_field());

	SECTION("velocity gauge"){
		perturbations::laser vector_potential({0.1, 0.0, 0.0}, 1.0_eV,perturbations::gauge::velocity);
		CHECK(vector_potential.has_uniform_vector_potential());
		CHECK(not vector_potential.has_uniform_electric_field());
		CHECK(vector_potential.uniform_vector_potential(0.0)[0] == 0.0);
		CHECK(vector_potential.uniform_vector_potential(0.0)[2] == 0.0);
	}

	SECTION("length gauge"){
		perturbations::laser E_field({1.0, 0.0, 0.0}, 1.0_eV);
		CHECK(E_field.has_uniform_electric_field());
		CHECK(not E_field.has_uniform_vector_potential());
		CHECK(E_field.uniform_electric_field(0.0)[0] == 0.0);
	}
}

#endif
#endif
