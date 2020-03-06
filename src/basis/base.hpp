/* -*- indent-tabs-mode: t -*- */

#ifndef BASIS__BASE
#define BASIS__BASE

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
#include <utils/partition.hpp>

namespace basis {

	/*
		This is a class that implements a very simple basis object. Useful for testing.
	*/
	
  class base {

  public:
		
		base(const long size, boost::mpi3::communicator & comm):
			part_(size, comm){
		}

		auto & part() {
			return part_;
		}

		auto & part() const {
			return part_;
		}
    
	protected:

		utils::partition part_;
		
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::base", "[basis::base]") {
  
  using namespace Catch::literals;
  using math::vec3d;
  
}
#endif

    
#endif
