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

		using element_type = Type;
		using basis_type = Basis;

		orbital_set(Basis const & basis, int const num_vectors, math::vector3<double> const & kpoint, boost::mpi3::cartesian_communicator<2> comm)
			:fields_(basis, num_vectors, comm),
       occupations_(fields_.local_set_size()),
			 kpoint_(kpoint){
		}
		
		orbital_set(Basis const & basis, int const num_vectors, boost::mpi3::cartesian_communicator<2> comm)
		:orbital_set(basis, num_vectors, math::vector3<double>{0.0, 0.0, 0.0}, comm){
		}
		
		orbital_set(basis::field_set<Basis, Type> && fields, math::array<double, 1> const & occs, math::vector3<double> const & kpoint)
		:fields_(std::move(fields)),
       occupations_(occs),
			 kpoint_(kpoint){
		}
		
    auto & fields() const {
      return fields_;
    }

    auto & fields() {
      return fields_;
    }

    auto & occupations() const {
      return occupations_;
    }

    auto & occupations() {
      return occupations_;
    }

		auto & kpoint() const {
			return kpoint_;
		}

		auto & set_part() const {
			return fields_.set_part();
		}

		auto local_set_size() const {
			return fields_.local_set_size();
		}
		
		auto set_size() const {
			return fields_.set_size();
		}
		
		auto & basis() const {
			return fields_.basis();
		}

		auto & basis() {
			return fields_.basis();
		}
		
		auto & matrix() const {
			return fields_.matrix();
		}

		auto & matrix() {
			return fields_.matrix();
		}
		
		auto cubic() const {
			return fields_.cubic();
		}

		auto cubic() {
			return fields_.cubic();
		}

		auto & full_comm() const {
			return fields_.full_comm();
		}

		auto & set_comm() const {
			return fields_.set_comm();
		}

	private:

    basis::field_set<Basis, Type> fields_;
    math::array<double, 1> occupations_;
		math::vector3<double> kpoint_;
		
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

	systems::box box = systems::box::orthorhombic(10.0_b, 4.0_b, 7.0_b).cutoff_energy(ecut);
  basis::real_space rs(box, basis_comm);

	states::orbital_set<basis::real_space, double> orb(rs, 12, cart_comm);

	CHECK(sizes(orb.fields().basis())[0] == 28);
	CHECK(sizes(orb.fields().basis())[1] == 11);
	CHECK(sizes(orb.fields().basis())[2] == 20);

	CHECK(orb.fields().local_set_size() == orb.local_set_size());
	CHECK(orb.fields().set_size() == orb.set_size());
	
	states::orbital_set<basis::real_space, double> orbk(rs, 12, {0.4, 0.22, -0.57}, cart_comm);

	CHECK(sizes(orbk.fields().basis())[0] == 28);
	CHECK(sizes(orbk.fields().basis())[1] == 11);
	CHECK(sizes(orbk.fields().basis())[2] == 20);

	CHECK(orbk.kpoint()[0] == 0.4_a);
	CHECK(orbk.kpoint()[1] == 0.22_a);
	CHECK(orbk.kpoint()[2] == -0.57_a);

	CHECK(orbk.fields().local_set_size() == orb.local_set_size());
	CHECK(orbk.fields().set_size() == orb.set_size());

	
}

#endif

#endif
