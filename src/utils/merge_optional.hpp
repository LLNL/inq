/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__MERGE_OPTIONAL
#define INQ__UTILS__MERGE_OPTIONAL

/*
 Copyright (C) 2019 Xavier Andrade

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
namespace utils {

template <class opt_type>
opt_type merge_optional(const opt_type & option1, const opt_type & option2){
	if(option2) return option2;
	if(option1) return option1;
	return opt_type{};
}

}
}
#endif

#ifdef INQ_UTILS_MERGE_OPTIONAL_UNIT_TEST
#undef INQ_UTILS_MERGE_OPTIONAL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
