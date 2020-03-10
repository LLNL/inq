/* -*- indent-tabs-mode: t -*- */

#ifndef BASIS_FIELD_SET
#define BASIS_FIELD_SET

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

#include <utils/partition.hpp>
#include <math/array.hpp>
#include <algorithm>
#ifdef HAVE_CUDA
#include <thrust/fill.h>
#endif

#include <mpi3/environment.hpp>
#include <mpi3/cartesian_communicator.hpp>
#include <utils/skeleton_wrapper.hpp>

namespace basis {
	
	template<class Basis, class type>
  class field_set {

  public:

		typedef Basis basis_type;
		typedef math::array<type, 2> internal_array_type;
		typedef type element_type;

		field_set(const basis_type & basis, const int num_vectors, boost::mpi3::cartesian_communicator<2> const & comm)
			:full_comm_(comm),
			 basis_comm_(comm.axis(1)),
			 set_comm_(comm.axis(0)),
			 set_part_(num_vectors, set_comm_),
			 matrix_({basis.part().local_size(), set_part_.local_size()}),
			 num_vectors_(num_vectors),
			 basis_(basis)
		{
			assert(basis_.part().comm_size() == basis_comm_.size());
    }

		field_set(const basis_type & basis, const int num_vectors, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance())
			:field_set(basis, num_vectors, boost::mpi3::cartesian_communicator<2>(comm))
		{
		}

		field_set(skeleton_wrapper<field_set<Basis, type>> const & skeleton)
			:field_set(skeleton.base.basis(), skeleton.base.set_size(), skeleton.base.basis_comm()){
		}

		auto skeleton(){
			return skeleton_wrapper<field_set<Basis, type>>(*this);
		}

		field_set(const field_set & coeff) = default;
		field_set(field_set && coeff) = default;
		field_set & operator=(const field_set & coeff) = default;
		field_set & operator=(field_set && coeff) = default;

		internal_array_type & matrix() {
			return matrix_;
		}

		internal_array_type const & matrix() const{
			return matrix_;
		}

		auto data() const {
			return matrix_.data();
		}

		auto data() {
			return matrix_.data();
		}

		auto num_elements() const {
			return matrix_.num_elements();
		}
		
		//set to a scalar value
		field_set & operator=(const type value) {

			//DATAOPERATIONS GPU::RUN FILL
			gpu::run(matrix_.num_elements(),
							 [lin = (element_type *) matrix_.data(), value] GPU_LAMBDA (auto ii){
								 lin[ii] = value;
							 });
			
			return *this;
		}
		
		const basis_type & basis() const {
			return basis_;
		}

		const int & set_size() const {
			return num_vectors_;
		}

		auto & set_part() const {
			return set_part_;
		}
		
		auto & set_comm() const {
			return set_comm_;
		}
				
		auto & basis_comm() const {
			return basis_comm_;
		}

		auto & full_comm() const {
			return full_comm_;
		}
		
		auto cubic() const {
			return matrix_.partitioned(basis_.cubic_dist(1).local_size()*basis_.cubic_dist(0).local_size()).partitioned(basis_.cubic_dist(0).local_size());
		}

		auto cubic() {
			return matrix_.partitioned(basis_.cubic_dist(1).local_size()*basis_.cubic_dist(0).local_size()).partitioned(basis_.cubic_dist(0).local_size());
		}

	private:

		mutable boost::mpi3::cartesian_communicator<2> full_comm_;
		mutable boost::mpi3::cartesian_communicator<1> basis_comm_;
		mutable boost::mpi3::cartesian_communicator<1> set_comm_;
		utils::partition set_part_;
		internal_array_type matrix_;
		int num_vectors_;
		basis_type basis_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

#include <mpi3/cartesian_communicator.hpp>

TEST_CASE("Class basis::field_set", "[basis::field_set]"){
  
  using namespace Catch::literals;
  using math::vec3d;
  
  double ecut = 40.0;

	auto comm = boost::mpi3::environment::get_world_instance();

	boost::mpi3::cartesian_communicator<2> cart_comm(comm);

	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

  ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 4.0, 0.0), vec3d(0.0, 0.0, 7.0));

  basis::real_space rs(cell, input::basis::cutoff_energy(ecut), basis_comm);

	basis::field_set<basis::real_space, double> ff(rs, 12, cart_comm);

	REQUIRE(sizes(rs)[0] == 28);
	REQUIRE(sizes(rs)[1] == 11);
	REQUIRE(sizes(rs)[2] == 20);

	//std::cout << ff.basis_comm().size() << " x " << ff.set_comm().size() << std::endl;
	//	std::cout << rs.part().comm_size() << std::endl;

	if(ff.basis_comm().size() == 1) REQUIRE(std::get<0>(sizes(ff.matrix())) == 6160);
	if(ff.basis_comm().size() == 2) REQUIRE(std::get<0>(sizes(ff.matrix())) == 6160/2);
	if(ff.set_comm().size() == 1) REQUIRE(std::get<1>(sizes(ff.matrix())) == 12);
	if(ff.set_comm().size() == 2) REQUIRE(std::get<1>(sizes(ff.matrix())) == 6);
	if(ff.set_comm().size() == 3) REQUIRE(std::get<1>(sizes(ff.matrix())) == 4);
	if(ff.set_comm().size() == 4) REQUIRE(std::get<1>(sizes(ff.matrix())) == 3);
	if(ff.set_comm().size() == 6) REQUIRE(std::get<1>(sizes(ff.matrix())) == 2);

	if(ff.basis_comm().size() == 1) REQUIRE(std::get<0>(sizes(ff.cubic())) == 28);
	if(ff.basis_comm().size() == 2) REQUIRE(std::get<0>(sizes(ff.cubic())) == 14);
	if(ff.basis_comm().size() == 4) REQUIRE(std::get<0>(sizes(ff.cubic())) == 7);
	REQUIRE(std::get<1>(sizes(ff.cubic())) == 11);
	REQUIRE(std::get<2>(sizes(ff.cubic())) == 20);
	if(ff.set_comm().size() == 1) REQUIRE(std::get<3>(sizes(ff.cubic())) == 12);
	if(ff.set_comm().size() == 2) REQUIRE(std::get<3>(sizes(ff.cubic())) == 6);

	ff = 12.2244;

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			REQUIRE(ff.matrix()[ii][jj] == 12.2244_a);
		}
	}
	
}

#endif

#endif
