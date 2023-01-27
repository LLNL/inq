/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__BASE
#define INQ__BASIS__BASE

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

#include <cassert>

#include <mpi3/environment.hpp>
#include <parallel/communicator.hpp>
#include <parallel/partition.hpp>

namespace inq {
namespace basis {

class base {
	
public:
	
	base(const long size, parallel::communicator & comm):
		comm_(comm),
		part_(size, comm){
	}
	
	auto & part() {
		return part_;
	}
	
	auto & part() const {
		return part_;
	}

	auto & comm() const {
		return comm_;
	}

	auto local_size() const {
		return part_.local_size();
	}

	protected:
	
	mutable parallel::communicator comm_;
	inq::parallel::partition part_;
		
};

}
}


#ifdef INQ_BASIS_BASE_UNIT_TEST
#undef INQ_BASIS_BASE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("class basis::base", "[basis::base]") {
  
	using namespace inq;
	using namespace Catch::literals;

}
#endif

#endif
