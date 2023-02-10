/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__RT
#define INPUT__RT

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <magnitude/time.hpp>
#include <utils/merge_optional.hpp>

#include <optional>
#include <cassert>


namespace inq {
namespace input {

class rt {

public:
	
	enum class electron_propagator { ETRS, CRANK_NICOLSON };
	
	rt(){
	}

	static auto dt(quantity<magnitude::time> dt) {
		rt solver;
		solver.dt_ = dt.in_atomic_units();
		return solver;
	}

	auto dt() const {
		return dt_.value_or(0.01);
	}

	auto static num_steps(double etol) {
		rt solver;
		solver.num_steps_ = etol;
		return solver;
	}
				
	auto num_steps() const {
		return num_steps_.value_or(100);
	}

	auto static etrs() {
		rt solver;
		solver.prop_ = electron_propagator::ETRS;
		return solver;
	}

	auto static crank_nicolson() {
		rt solver;
		solver.prop_ = electron_propagator::CRANK_NICOLSON;
		return solver;
	}
	
	auto propagator() const {
		return prop_.value_or(electron_propagator::ETRS);
	}
	
	friend auto operator|(const rt & solver1, const rt & solver2){
		using utils::merge_optional;

		rt rsolver;
		rsolver.dt_	= merge_optional(solver1.dt_, solver2.dt_);
		rsolver.num_steps_	= merge_optional(solver1.num_steps_, solver2.num_steps_);
		rsolver.prop_	= merge_optional(solver1.prop_, solver2.prop_);		
		return rsolver;
	}
    
private:

	std::optional<double> dt_;
	std::optional<int> num_steps_;
	std::optional<electron_propagator> prop_;
		
};
    
}
}
#endif

#ifdef INQ_INPUT_RT_UNIT_TEST
#undef INQ_INPUT_RT_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

  using namespace inq;
  using namespace inq::magnitude;	
  using namespace Catch::literals;

	SECTION("Defaults"){

    input::rt solver;

    CHECK(solver.dt() == 0.01_a);
    CHECK(solver.num_steps() == 100);
    CHECK(solver.propagator() == input::rt::electron_propagator::ETRS);		
    
  }

  SECTION("Composition"){

    auto solver = input::rt::num_steps(1000) | input::rt::dt(0.05_atomictime) | input::rt::crank_nicolson();
    
    CHECK(solver.num_steps() == 1000);
    CHECK(solver.dt() == 0.05_a);
		CHECK(solver.propagator() == input::rt::electron_propagator::CRANK_NICOLSON);
  }

}
#endif
