/* -*- indent-tabs-mode: t -*- */

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

#ifndef INQ__INQ
#define INQ__INQ

#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <input/parse_xyz.hpp>
#include <operations/io.hpp>
#include <real_time/propagate.hpp>
#include <utils/match.hpp>
#include <perturbations/kick.hpp>
#include <perturbations/laser.hpp>
#include <observables/spectrum.hpp>

#ifdef INQ_INQ_INQ_UNIT_TEST
#undef INQ_INQ_INQ_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("inq::inq", "[inq::inq]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif
