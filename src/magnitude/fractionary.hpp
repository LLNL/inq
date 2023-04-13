/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__FRACTIONARY
#define INQ__MAGNITUDE__FRACTIONARY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/quantity.hpp>

namespace inq {
namespace magnitude {
    
class fractionary {
    
};
  
auto operator "" _fractionary(long double val){
  return inq::quantity<fractionary>::from_atomic_units(val);
}
  
auto operator "" _frac(long double val){
  return inq::quantity<fractionary>::from_atomic_units(val);
}

auto operator "" _relative(long double val){
  return inq::quantity<fractionary>::from_atomic_units(val);
}

auto operator "" _rel(long double val){
  return inq::quantity<fractionary>::from_atomic_units(val);
}  

auto operator "" _crystal(long double val){
  return inq::quantity<fractionary>::from_atomic_units(val);
}

auto operator "" _crys(long double val){
  return inq::quantity<fractionary>::from_atomic_units(val);
}  

GPU_FUNCTION auto in_atomic_units(vector3<quantity<magnitude::fractionary>> const & quant) {
	return vector3<double, contravariant>{quant[0].in_atomic_units(), quant[1].in_atomic_units(), quant[2].in_atomic_units()};
}


}
}
#endif

#ifdef INQ_MAGNITUDE_FRACTIONARY_UNIT_TEST
#undef INQ_MAGNITUDE_FRACTIONARY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using namespace magnitude;
  
}
#endif

