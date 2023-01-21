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
		using kpoint_type = math::vector3<double, math::covariant>;

		orbital_set(Basis const & basis, int const num_vectors, kpoint_type const & kpoint, int spin_index, parallel::cartesian_communicator<2> comm)
			:fields_(basis, num_vectors, comm),
			 kpoint_(kpoint),
			 spin_index_(spin_index){
		}
		
		orbital_set(orbital_set && oldset, parallel::cartesian_communicator<2> new_comm)
		:fields_(std::move(oldset.fields_), new_comm),
			 kpoint_(oldset.kpoint()),
			 spin_index_(oldset.spin_index()){
		}
		
		template <class any_type>
		orbital_set(inq::utils::skeleton_wrapper<orbital_set<Basis, any_type>> const & skeleton)
			:orbital_set(skeleton.base.basis(), skeleton.base.set_size(), skeleton.base.kpoint(), skeleton.base.spin_index(), skeleton.base.full_comm()){
		}
		
		auto skeleton() const {
			return inq::utils::skeleton_wrapper<orbital_set<Basis, Type>>(*this);
		}

		template <class OtherType>
		static auto reciprocal(inq::utils::skeleton_wrapper<orbital_set<typename basis_type::reciprocal_space, OtherType>> const & skeleton){
			return orbital_set<basis_type, element_type>(skeleton.base.basis().reciprocal(), skeleton.base.set_size(),  skeleton.base.kpoint(), skeleton.base.spin_index(), skeleton.base.full_comm());
		}

		template <typename ScalarType>
		auto fill(ScalarType const & scalar){
			fields_.fill(scalar);			
		}
		
		auto & kpoint() const {
			return kpoint_;
		}

		auto & spin_index() const {
			assert(spin_index_ >= 0 and spin_index_ < 2);
			return spin_index_;
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

		auto hypercubic() const {
			return fields_.hypercubic();
		}

		auto hypercubic() {
			return fields_.hypercubic();
		}

		auto & full_comm() const {
			return fields_.full_comm();
		}

		auto & set_comm() const {
			return fields_.set_comm();
		}

		auto par_set_begin() const {
			return fields_.par_set_begin();
		}

		auto par_set_end() const {
			return fields_.par_set_end();
		}
		
	private:

    basis::field_set<Basis, Type> fields_;
		kpoint_type kpoint_;
		int spin_index_;
		
  };

}
}

#ifdef INQ_STATES_ORBITAL_SET_UNIT_TEST
#undef INQ_STATES_ORBITAL_SET_UNIT_TEST

#include <basis/real_space.hpp>

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>

#include <parallel/communicator.hpp>

TEST_CASE("Class states::orbital_set", "[states::orbital_set]"){
  
	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
  using math::vector3;
  
  auto ecut = 40.0_Ha;

	auto comm = boost::mpi3::environment::get_world_instance();

	parallel::cartesian_communicator<2> cart_comm(comm, {});

	auto set_comm = basis::set_subcomm(cart_comm);
	auto basis_comm = basis::basis_subcomm(cart_comm);	

	systems::box box = systems::box::orthorhombic(10.0_b, 4.0_b, 7.0_b).cutoff_energy(ecut);
  basis::real_space rs(box, basis_comm);

	states::orbital_set<basis::real_space, double> orb(rs, 12, math::vector3<double, math::covariant>{0.0, 0.0, 0.0}, 0, cart_comm);

	CHECK(sizes(orb.basis())[0] == 28);
	CHECK(sizes(orb.basis())[1] == 11);
	CHECK(sizes(orb.basis())[2] == 20);

	CHECK(orb.local_set_size() == orb.local_set_size());
	CHECK(orb.set_size() == orb.set_size());
	
	states::orbital_set<basis::real_space, double> orbk(rs, 12, {0.4, 0.22, -0.57}, 0, cart_comm);

	CHECK(sizes(orbk.basis())[0] == 28);
	CHECK(sizes(orbk.basis())[1] == 11);
	CHECK(sizes(orbk.basis())[2] == 20);

	CHECK(orbk.kpoint()[0] == 0.4_a);
	CHECK(orbk.kpoint()[1] == 0.22_a);
	CHECK(orbk.kpoint()[2] == -0.57_a);

	CHECK(orbk.local_set_size() == orb.local_set_size());
	CHECK(orbk.set_size() == orb.set_size());

	states::orbital_set<basis::real_space, double> orb_copy(orbk.skeleton());

	CHECK(sizes(orb_copy.basis()) == sizes(orbk.basis()));
	CHECK(orb_copy.kpoint() == orbk.kpoint());
	CHECK(orb_copy.local_set_size() == orbk.local_set_size());
	CHECK(orb_copy.set_size() == orbk.set_size());
	
}

#endif

#endif
