/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__FIELD_SET
#define INQ__BASIS__FIELD_SET

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

#include <basis/real_space.hpp>
#include <utils/partition.hpp>
#include <math/array.hpp>
#include <algorithm>

#include <mpi3/environment.hpp>
#include <mpi3/cartesian_communicator.hpp>
#include <utils/skeleton_wrapper.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace basis {
	
	template<class Basis, class type>
  class field_set {

  public:

		typedef Basis basis_type;
		typedef math::array<type, 2> internal_array_type;
		typedef type element_type;

		field_set(const basis_type & basis, const int num_vectors, boost::mpi3::cartesian_communicator<2> const & comm)
			:full_comm_(comm),
			 set_comm_(comm.axis(0)),
			 set_part_(num_vectors, set_comm_),
			 matrix_({basis.part().local_size(), set_part_.local_size()}),
			 num_vectors_(num_vectors),
			 basis_(basis)
		{
			prefetch();
			assert(basis_.part().comm_size() == comm.axis(1).size());
			assert(local_set_size() > 0);
    }

		//when no communicator is given, use the basis communicator
		field_set(const basis_type & basis, const int num_vectors)			
			:field_set(basis, num_vectors, boost::mpi3::cartesian_communicator<2>(basis.comm(), {1}))
		{
		}

		template <class any_type>
		field_set(inq::utils::skeleton_wrapper<field_set<Basis, any_type>> const & skeleton)
			:field_set(skeleton.base.basis(), skeleton.base.set_size(), skeleton.base.full_comm()){
		}
		
		// Avoid the default copy constructor since the multi copy constructor is slow
		//		field_set(const field_set & coeff) = default;
		field_set(field_set const & other)
			:field_set(other.skeleton()){
			matrix_ = other.matrix_;
		}

		field_set(field_set && coeff) = default;
		field_set & operator=(const field_set & coeff) = default;
		field_set & operator=(field_set && coeff) = default;

		auto skeleton() const {
			return inq::utils::skeleton_wrapper<field_set<Basis, type>>(*this);
		}
	
		internal_array_type & matrix() {
			return matrix_;
		}

		internal_array_type const & matrix() const{
			return matrix_;
		}

		auto data() const {
			return raw_pointer_cast(matrix_.data());
		}

		auto data() {
			return raw_pointer_cast(matrix_.data_elements());
		}

		auto num_elements() const {
			return matrix_.num_elements();
		}
		
		//set to a scalar value
		field_set & operator=(const type value) {

			CALI_CXX_MARK_SCOPE("field_set=scalar");

			//DATAOPERATIONS GPU::RUN
			gpu::run(matrix_.num_elements(),
							 [lin = (element_type *) matrix_.data_elements(), value] GPU_LAMBDA (auto ii){
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

		auto local_set_size() const {
			return set_part_.local_size();
		}
		
		auto & set_part() const {
			return set_part_;
		}
		
		auto & set_comm() const {
			return set_comm_;
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

		void prefetch() const {
			math::prefetch(matrix_);
		}
		
	private:

		mutable boost::mpi3::cartesian_communicator<2> full_comm_;
		mutable boost::mpi3::cartesian_communicator<1> set_comm_;
		inq::utils::partition set_part_;
		internal_array_type matrix_;
		int num_vectors_;
		basis_type basis_;

  };

auto complex_field(field_set<basis::real_space, double> const & field) {
	field_set<basis::real_space, inq::complex> cfield(field.skeleton());
	cfield.matrix() = field.matrix();
	return cfield;
}

field_set<basis::real_space, double> real_field(field_set<basis::real_space, inq::complex> const & field) {
	
	field_set<basis::real_space, double> rfield(field.skeleton());
	
	// Multi should be able to do this, but it causes a lot of compilation troubles
	//			rfield.matrix() = boost::multi::blas::real(matrix());
	
	//DATAOPERATIONS GPU::RUN 1D
	gpu::run(field.set_part().local_size(), field.basis().part().local_size(),
					 [rp = begin(rfield.matrix()), cp = begin(field.matrix())] GPU_LAMBDA (auto ist, auto ii){
						 rp[ii][ist] = inq::real(cp[ii][ist]);
					 });
	
	return rfield;
}

}
}

#ifdef INQ_BASIS_FIELD_SET_UNIT_TEST
#undef INQ_BASIS_FIELD_SET_UNIT_TEST

#include <basis/real_space.hpp>

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

#include <mpi3/cartesian_communicator.hpp>

TEST_CASE("Class basis::field_set", "[basis::field_set]"){
  
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

	basis::field_set<basis::real_space, double> ff(rs, 12, cart_comm);

	CHECK(sizes(rs)[0] == 28);
	CHECK(sizes(rs)[1] == 11);
	CHECK(sizes(rs)[2] == 20);

	//std::cout << ff.basis().comm().size() << " x " << ff.set_comm().size() << std::endl;
	//	std::cout << rs.part().comm_size() << std::endl;

	if(ff.basis().comm().size() == 1) CHECK(std::get<0>(sizes(ff.matrix())) == 6160);
	if(ff.basis().comm().size() == 2) CHECK(std::get<0>(sizes(ff.matrix())) == 6160/2);
	if(ff.set_comm().size() == 1) CHECK(std::get<1>(sizes(ff.matrix())) == 12);
	if(ff.set_comm().size() == 2) CHECK(std::get<1>(sizes(ff.matrix())) == 6);
	if(ff.set_comm().size() == 3) CHECK(std::get<1>(sizes(ff.matrix())) == 4);
	if(ff.set_comm().size() == 4) CHECK(std::get<1>(sizes(ff.matrix())) == 3);
	if(ff.set_comm().size() == 6) CHECK(std::get<1>(sizes(ff.matrix())) == 2);

	if(ff.basis().comm().size() == 1) CHECK(std::get<0>(sizes(ff.cubic())) == 28);
	if(ff.basis().comm().size() == 2) CHECK(std::get<0>(sizes(ff.cubic())) == 14);
	if(ff.basis().comm().size() == 4) CHECK(std::get<0>(sizes(ff.cubic())) == 7);
	CHECK(std::get<1>(sizes(ff.cubic())) == 11);
	CHECK(std::get<2>(sizes(ff.cubic())) == 20);
	if(ff.set_comm().size() == 1) CHECK(std::get<3>(sizes(ff.cubic())) == 12);
	if(ff.set_comm().size() == 2) CHECK(std::get<3>(sizes(ff.cubic())) == 6);

	ff = 12.2244;

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(ff.matrix()[ii][jj] == 12.2244_a);
		}
	}

	auto zff = complex_field(ff);
	
	static_assert(std::is_same<decltype(zff), basis::field_set<basis::real_space, complex>>::value, "complex() should return a complex field");
		
	CHECK(std::get<1>(sizes(zff.cubic())) == 11);
	CHECK(std::get<2>(sizes(zff.cubic())) == 20);

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(real(zff.matrix()[ii][jj]) == 12.2244_a);
			CHECK(imag(zff.matrix()[ii][jj]) == 0.0_a);
		}
	}

	auto dff = real_field(zff);

	static_assert(std::is_same<decltype(dff), basis::field_set<basis::real_space, double>>::value, "real() should return a double field");

	CHECK(std::get<1>(sizes(dff.cubic())) == 11);
	CHECK(std::get<2>(sizes(dff.cubic())) == 20);

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(dff.matrix()[ii][jj] == 12.2244_a);
		}
	}

	dff = 0.0;

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(dff.matrix()[ii][jj] == 0.0_a);
		}
	}
	
	zff = 0.0;
	
	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(real(zff.matrix()[ii][jj]) == 0.0_a);
			CHECK(imag(zff.matrix()[ii][jj]) == 0.0_a);
		}
	}
	
}

#endif

#endif
