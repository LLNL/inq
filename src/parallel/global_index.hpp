/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__GLOBAL_INDEX
#define INQ__PARALLEL__GLOBAL_INDEX

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


#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

class global_index {

 public:

	constexpr explicit global_index(long val):
		value_(val){
		}

	constexpr auto & value() const {
		return value_;
	}
	
 private:
	long value_; 
};

}
}

#ifdef INQ_PARALLEL_GLOBAL_INDEX_UNIT_TEST
#undef INQ_PARALLEL_GLOBAL_INDEX_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class parallel::global_index", "[parallel::global_index]") {
  
	using namespace inq;
	using namespace Catch::literals;

	parallel::global_index gi(10);

	CHECK(gi.value() == 10);
	
}
#endif

    
#endif
