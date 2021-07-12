/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__STATES__ORBITAL_SET
#define INQ__STATES__ORBITAL_SET

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

#include <basis/field_set.hpp>

namespace inq {
namespace states {
	
	template<class Basis, class Type>
  class orbital_set {

  public:

		orbital_set(Basis const & basis, int const num_vectors, boost::mpi3::cartesian_communicator<2> comm)
			:fields_(basis, num_vectors, comm){
		}

    auto & fields() const {
      return fields_;
    }

    auto & fields() {
      return fields_;
    }
    
	private:

    basis::field_set<Basis, Type> fields_;
    
  };

}
}

#ifdef INQ_STATES_ORBITAL_SET_UNIT_TEST
#undef INQ_STATES_ORBITAL_SET_UNIT_TEST

#include <basis/real_space.hpp>

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

#include <mpi3/cartesian_communicator.hpp>

TEST_CASE("Class states::orbital_set", "[states::orbital_set]"){
  
	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
  using math::vector3;
  
  auto ecut = 40.0_Ha;

	auto comm = boost::mpi3::environment::get_world_instance();

	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {});

	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

  ions::UnitCell cell(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 4.0, 0.0), vector3<double>(0.0, 0.0, 7.0));

  basis::real_space rs(cell, input::basis::cutoff_energy(ecut), basis_comm);

	states::orbital_set<basis::real_space, double> orb(rs, 12, cart_comm);

	CHECK(sizes(orb.fields().basis())[0] == 28);
	CHECK(sizes(orb.fields().basis())[1] == 11);
	CHECK(sizes(orb.fields().basis())[2] == 20);
	
}

#endif

#endif
