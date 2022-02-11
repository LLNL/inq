/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__FRACTIONARY
#define INQ__MAGNITUDE__FRACTIONARY

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


}
}

#ifdef INQ_MAGNITUDE_FRACTIONARY_UNIT_TEST
#undef INQ_MAGNITUDE_FRACTIONARY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("magnitude::fractionary", "[magnitude::fractionary]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace magnitude;

  
}

#endif

#endif

