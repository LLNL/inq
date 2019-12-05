/* -*- indent-tabs-mode: t -*- */

#ifndef UTILS__DISTRIBUTION
#define UTILS__DISTRIBUTION

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
#include <array>

namespace utils {

  template <class comm_type>
  class distribution {

  public:

		distribution(const comm_type & comm):
      comm_(comm){
		}
    
	protected:
    
    comm_type comm_;
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("class utils::distribution", "[utils::distribution]") {
  
  using namespace Catch::literals;
  using math::d3vector;

  utils::distribution<boost::mpi3::communicator> dist(boost::mpi3::environment::get_world_instance());
  
}
#endif

    
#endif
