/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__STATES__INDEX
#define INQ__STATES__INDEX

/*
 Copyright (C) 2023 Xavier Andrade

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

#include <math/vector3.hpp>

#include <unordered_map>

namespace inq {
namespace states {

struct key {
  vector3<double, covariant> kpoint;
  int spin_index;

  bool operator==(key const & other) const {
    return kpoint == other.kpoint and spin_index == other.spin_index;
  }
};

struct hash {
  auto operator()(key const & kk) const {
    auto h0 = std::hash<int>{}(kk.spin_index);
    h0 ^= (std::hash<double>{}(kk.kpoint[0]) << 1);
    h0 ^= (std::hash<double>{}(kk.kpoint[1]) << 1);
    h0 ^= (std::hash<double>{}(kk.kpoint[2]) << 1);
    return h0;
  }
};

using index = std::unordered_map<key, int, hash>;

}
}
#endif

#ifdef INQ_STATES_INDEX_UNIT_TEST
#undef INQ_STATES_INDEX_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <mpi3/environment.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace Catch::literals;
	using namespace inq;

  states::key key1{vector3<double, covariant>(1.0, 2.0, 3.0), 4};
  states::key key2{vector3<double, covariant>(0.56, -2.0, 12.4), 0};
  states::key key3{vector3<double, covariant>(2.33, 0.44, 3.9), 0};

  CHECK(key1 == key1);
  
  CHECK(states::hash{}(key1) != states::hash{}(key2));
  CHECK(states::hash{}(key1) != states::hash{}(key3));
  
  states::index idx;

  idx[key1] = 1;
  idx[key2] = 2;
  idx[key3] = 3;

  CHECK(idx[key1] == 1);
  CHECK(idx[key2] == 2);
  CHECK(idx[key3] == 3);

  CHECK(idx.size() == 3);
    
}
#endif
