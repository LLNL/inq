/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARALLEL__COMMUNICATOR
#define INQ__PARALLEL__COMMUNICATOR

/*
 Copyright (C) 2022 Xavier Andrade

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
#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/environment.hpp>

#include <cassert>
#include <array>

namespace inq{
namespace parallel {

using communicator = boost::mpi3::communicator;

template<boost::mpi3::dimensionality_type D = boost::mpi3::dynamic_extent>
using cartesian_communicator = boost::mpi3::cartesian_communicator<D>;

}
}

#ifdef INQ_PARALLEL_COMMUNICATOR_UNIT_TEST
#undef INQ_PARALLEL_COMMUNICATOR_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("class parallel::communicator", "[parallel::communicator]") {
  
}
#endif

    
#endif
