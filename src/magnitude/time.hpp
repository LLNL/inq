/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MAGNITUDE__TIME
#define INQ__MAGNITUDE__TIME

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
#include <magnitude/energy.hpp>

namespace inq {
namespace magnitude {

class time {
	
};

auto operator "" _atomictime(long double val){
	return inq::quantity<time>::from_atomic_units(val);
}

auto operator "" _attosecond(long double val){
	return inq::quantity<time>::from_atomic_units(0.0413413733352975*val);
}

auto operator "" _as(long double val){
	return val*1.0_attosecond;
}

auto operator "" _femtosecond(long double val){
	return val*1000.0_attosecond;
}

auto operator "" _fs(long double val){
	return val*1.0_femtosecond;
}

auto operator "" _picosecond(long double val){
	return val*1000.0_femtosecond;
}

auto operator "" _ps(long double val){
	return val*1.0_picosecond;
}

auto operator "" _nanosecond(long double val){
	return val*1000.0_picosecond;
}

auto operator "" _ns(long double val){
	return val*1.0_nanosecond;
}

auto operator/(double num, quantity<energy> den){
  return quantity<time>::from_atomic_units(num/den.in_atomic_units());
}

}
}
#endif

#ifdef INQ_MAGNITUDE_TIME_UNIT_TEST
#undef INQ_MAGNITUDE_TIME_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("magnitude::time", "[magnitude::time]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace magnitude;

  {
    auto ti = 100.0_atomictime;
    CHECK(ti.in_atomic_units() == 100.0_a);
  }

  {
    auto ti = 0.9562_attosecond;
    CHECK(ti.in_atomic_units() == 0.0395306211832115_a);
  }
  
  {
    auto ti = 0.9562_as;
    CHECK(ti.in_atomic_units() == 0.0395306211832115_a);
  }

  {
    auto ti = 43.27_femtosecond;
    CHECK(ti.in_atomic_units() == 1788.84122421832_a);
  }
   
  {
    auto ti = 43.27_fs;
    CHECK(ti.in_atomic_units() == 1788.84122421832_a);
  }
  
  {
    auto ti = 17.77_picosecond;
    CHECK(ti.in_atomic_units() == 734636.204168237_a);
  }
  
  {
    auto ti = 17.77_ps;
    CHECK(ti.in_atomic_units() == 734636.204168237_a);
  }

  {
    auto ti = 0.03974_nanosecond;
    CHECK(ti.in_atomic_units() == 1642906.17634472_a);
  }

  {
    auto ti = 0.03974_ns;
    CHECK(ti.in_atomic_units() == 1642906.17634472_a);
  }
  
  {
    auto ti = 1.0/1.0_Ha;
    CHECK(ti.in_atomic_units() == 1.0_a);
  }

  {
    auto ti = 1.0/100.0_eV;
    CHECK(ti.in_atomic_units() == 0.272113862460642_a);
  }
  
}
#endif

