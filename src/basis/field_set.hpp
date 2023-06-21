/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__FIELD_SET
#define INQ__BASIS__FIELD_SET

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <basis/real_space.hpp>
#include <parallel/partition.hpp>
#include <parallel/get_remote_points.hpp>
#include <gpu/array.hpp>
#include <algorithm>

#include <mpi3/environment.hpp>
#include <parallel/communicator.hpp>
#include <utils/skeleton_wrapper.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <gpu/copy.hpp>

namespace inq {
namespace basis {

auto set_subcomm(parallel::cartesian_communicator<2> & comm){
	return parallel::cartesian_communicator<1>(comm.axis(1));
}
auto basis_subcomm(parallel::cartesian_communicator<2> & comm){
	return parallel::cartesian_communicator<1>(comm.axis(0));
}

template<class BasisType, class ElementType, class PartitionType = inq::parallel::partition>
class field_set {

public:

	typedef BasisType basis_type;
	typedef gpu::array<ElementType, 2> internal_array_type;
	typedef ElementType element_type;

	template <typename BType, typename EType>
	using template_type = field_set<BType, EType>;
	
	field_set(const basis_type & basis, PartitionType part, parallel::cartesian_communicator<2> comm)
		:full_comm_(std::move(comm)),
		 set_comm_(basis::set_subcomm(full_comm_)),
		 set_part_(std::move(part)),
		 matrix_({basis.part().local_size(), set_part_.local_size()}),
		 num_vectors_(set_part_.size()),
		 basis_(basis)
	{
		prefetch();
		assert(part.comm_size() == basis::set_subcomm(full_comm_).size());
		assert(basis_.part().comm_size() == basis::basis_subcomm(full_comm_).size());
		assert(local_set_size() > 0);
	}

	field_set(const basis_type & basis, const int num_vectors, parallel::cartesian_communicator<2> comm)
		:field_set(basis, parallel::partition(num_vectors, set_subcomm(comm)), comm){
	}

	//when no communicator is given, use the basis communicator
	field_set(const basis_type & basis, const int num_vectors)			
		:field_set(basis, num_vectors, parallel::cartesian_communicator<2>(basis.comm(), {basis.comm().size(), 1}))
	{
	}

	template <class AnyType, class AnyPartType>
	field_set(inq::utils::skeleton_wrapper<field_set<BasisType, AnyType, AnyPartType>> const & skeleton)
		:field_set(skeleton.base.basis(), skeleton.base.set_size(), skeleton.base.full_comm()){
	}
		
	// Avoid the default copy constructor since the multi copy constructor is slow
	//		field_set(const field_set & coeff) = default;
	field_set(field_set const & other)
		:field_set(other.skeleton()){
		matrix_ = other.matrix_;
	}

	field_set(field_set && coeff) = default;

	field_set(field_set && oldset, parallel::cartesian_communicator<2> new_comm):
		field_set(BasisType{BasisType{oldset.basis()}, basis_subcomm(new_comm)}, oldset.set_size(), new_comm)
	{
		gpu::array<int, 1> rem_points(basis().local_size());
		gpu::array<int, 1> rem_states(local_set_size());
		for(long ip = 0; ip < basis().local_size(); ip++) rem_points[ip] = basis().part().local_to_global(ip).value();
		for(long ist = 0; ist < local_set_size(); ist++) rem_states[ist] = set_part().local_to_global(ist).value();
		matrix_ = parallel::get_remote_points(oldset, rem_points, rem_states);
	}
				
	field_set & operator=(field_set const& other){
		full_comm_ = other.full_comm_;
		set_comm_  = other.set_comm_;
		set_part_  = other.set_part_;
		matrix_    = other.matrix_;
		num_vectors_ = other.num_vectors_;
		basis_     = other.basis_;
		return *this;
	}
	field_set & operator=(field_set && coeff) = default;

	auto skeleton() const {
		return inq::utils::skeleton_wrapper<field_set<BasisType, ElementType, PartitionType>>(*this);
	}

	template <class OtherType, class AnyPartType>
	static auto reciprocal(inq::utils::skeleton_wrapper<field_set<basis_type, OtherType, AnyPartType>> const & skeleton){
		return field_set<typename basis_type::reciprocal_space, element_type, PartitionType>(skeleton.base.basis().reciprocal(), skeleton.base.set_size(), skeleton.base.full_comm());
	}
	
	internal_array_type & matrix() {
		return matrix_;
	}

	internal_array_type const & matrix() const{
		return matrix_;
	}

	auto data() const {
		return raw_pointer_cast(matrix_.data_elements());
	}

	auto data() {
		return raw_pointer_cast(matrix_.data_elements());
	}

	auto num_elements() const {
		return matrix_.num_elements();
	}

	template <typename ScalarType>
	void fill(ScalarType const & scalar) {
		CALI_CXX_MARK_SCOPE("fill(field_set)");

		gpu::run(matrix_.num_elements(), [lin = raw_pointer_cast(matrix_.data_elements()), scalar] GPU_LAMBDA (auto ii){
			lin[ii] = scalar;
		});
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

	auto hypercubic() const {
		return matrix_.partitioned(basis_.cubic_part(1).local_size()*basis_.cubic_part(0).local_size()).partitioned(basis_.cubic_part(0).local_size());
	}

	auto hypercubic() {
		return matrix_.partitioned(basis_.cubic_part(1).local_size()*basis_.cubic_part(0).local_size()).partitioned(basis_.cubic_part(0).local_size());
	}

	void prefetch() const {
		gpu::prefetch(matrix_);
	}

	class parallel_set_iterator {
			
		internal_array_type matrix_;
		int istep_;
		mutable parallel::cartesian_communicator<1> set_comm_;
		PartitionType set_part_;
		
	public:
			
		parallel_set_iterator(long basis_local_size, PartitionType set_part, parallel::cartesian_communicator<1> set_comm, internal_array_type const & data):
			matrix_({basis_local_size, set_part.max_local_size()}),
			istep_(0),
			set_comm_(std::move(set_comm)),
			set_part_(std::move(set_part)){

			CALI_CXX_MARK_SCOPE("field_set_iterator_constructor");
 
			gpu::copy(basis_local_size, set_part.local_size(), data, matrix_);
			set_comm_.nccl_init();
		};
		
		void operator++(){

			CALI_CXX_MARK_SCOPE("field_set_iterator++");
				
			auto mpi_type = boost::mpi3::detail::basic_datatype<element_type>();
				
			auto next_proc = set_comm_.rank() + 1;
			if(next_proc == set_comm_.size()) next_proc = 0;
			auto prev_proc = set_comm_.rank() - 1;
			if(prev_proc == -1) prev_proc = set_comm_.size() - 1;

			if(istep_ < set_comm_.size() - 1) {  //there is no need to copy for the last step

#ifdef ENABLE_NCCL
				ncclGroupStart();
				auto copy = matrix_;
				ncclRecv(raw_pointer_cast(matrix_.data_elements()), matrix_.num_elements()*sizeof(ElementType)/sizeof(double), ncclDouble, next_proc, &set_comm_.nccl_comm(), 0);
				ncclSend(raw_pointer_cast(copy.data_elements()), matrix_.num_elements()*sizeof(ElementType)/sizeof(double), ncclDouble, prev_proc, &set_comm_.nccl_comm(), 0);
				ncclGroupEnd();
				gpu::sync();
#else
				MPI_Sendrecv_replace(raw_pointer_cast(matrix_.data_elements()), matrix_.num_elements(), mpi_type, prev_proc, istep_, next_proc, istep_, set_comm_.get(), MPI_STATUS_IGNORE);
#endif
			}
				
			istep_++;
		}

		bool operator!=(int it_istep){
			return istep_ != it_istep;
		}

		auto matrix() const {
			return matrix_(boost::multi::ALL, {0, set_part_.local_size(set_ipart())});
		}

		auto set_ipart() const {
			auto ip = istep_ + set_comm_.rank();
			if(ip >= set_comm_.size()) ip -= set_comm_.size();
			return ip;
		}
			
	};

	auto par_set_begin() const {
		return parallel_set_iterator(basis().local_size(), set_part_, set_comm_, matrix());
	}

	auto par_set_end() const {
		return set_comm_.size();
	}

	template <typename CommunicatorType, typename OpType = std::plus<>>
	void all_reduce(CommunicatorType & comm, OpType op = OpType{}){
		if(comm.size() < 2) return;
		comm.all_reduce_n(raw_pointer_cast(matrix().data_elements()), matrix().num_elements(), op);
	}
		
private:

	mutable parallel::cartesian_communicator<2> full_comm_;
	mutable parallel::cartesian_communicator<1> set_comm_;
	PartitionType set_part_;
	internal_array_type matrix_;
	int num_vectors_;
	basis_type basis_;

};

field_set<basis::real_space, inq::complex> complex_field(field_set<basis::real_space, double> const & field) {
	field_set<basis::real_space, inq::complex> cfield(field.skeleton());

	gpu::run(field.set_part().local_size(), field.basis().part().local_size(),
					 [fie = begin(field.matrix()), cfie = begin(cfield.matrix())] GPU_LAMBDA (auto ist, auto ii){
						 cfie[ii][ist] = complex(fie[ii][ist], 0.0);
					 });
	
	return cfield;
}

template <class VectorSpace>
field_set<basis::real_space, vector3<complex, VectorSpace>> complex_field(field_set<basis::real_space, vector3<double, VectorSpace>> const & field) {
	field_set<basis::real_space, vector3<complex, VectorSpace>> cfield(field.skeleton());

	gpu::run(field.set_part().local_size(), field.basis().part().local_size(),
					 [fie = begin(field.matrix()), cfie = begin(cfield.matrix())] GPU_LAMBDA (auto ist, auto ii){
						 cfie[ii][ist] = complex{1.0, 0.0}*fie[ii][ist];
					 });
	
	return cfield;
}

field_set<basis::real_space, double> real_field(field_set<basis::real_space, inq::complex> const & field) {
	
	field_set<basis::real_space, double> rfield(field.skeleton());
	
	// Multi should be able to do this, but it causes a lot of compilation troubles
	//			rfield.matrix() = boost::multi::blas::real(matrix());
	
	gpu::run(field.set_part().local_size(), field.basis().part().local_size(),
					 [rp = begin(rfield.matrix()), cp = begin(field.matrix())] GPU_LAMBDA (auto ist, auto ii){
						 rp[ii][ist] = inq::real(cp[ii][ist]);
					 });
	
	return rfield;
}

template <class VectorSpace>
field_set<basis::real_space, vector3<double, VectorSpace>> real_field(field_set<basis::real_space, vector3<complex, VectorSpace>> const & field) {
	
	field_set<basis::real_space, vector3<double, VectorSpace>> rfield(field.skeleton());
	
	gpu::run(field.set_part().local_size(), field.basis().part().local_size(),
					 [rp = begin(rfield.matrix()), cp = begin(field.matrix())] GPU_LAMBDA (auto ist, auto ii){
						 rp[ii][ist] = real(cp[ii][ist]);
					 });
	
	return rfield;
}

}
}
#endif

#ifdef INQ_BASIS_FIELD_SET_UNIT_TEST
#undef INQ_BASIS_FIELD_SET_UNIT_TEST

#include <basis/real_space.hpp>

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>

#include <parallel/communicator.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){
  
	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	parallel::cartesian_communicator<2> cart_comm(comm, {});

	auto set_comm = basis::set_subcomm(cart_comm);
	auto basis_comm = basis::basis_subcomm(cart_comm);	

	systems::box box = systems::box::orthorhombic(10.0_b, 4.0_b, 7.0_b);
	basis::real_space rs(box, /*spacing = */ 0.35124074, basis_comm);
	
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

	if(ff.basis().comm().size() == 1) CHECK(std::get<0>(sizes(ff.hypercubic())) == 28);
	if(ff.basis().comm().size() == 2) CHECK(std::get<0>(sizes(ff.hypercubic())) == 14);
	if(ff.basis().comm().size() == 4) CHECK(std::get<0>(sizes(ff.hypercubic())) == 7);
	CHECK(std::get<1>(sizes(ff.hypercubic())) == 11);
	CHECK(std::get<2>(sizes(ff.hypercubic())) == 20);
	if(ff.set_comm().size() == 1) CHECK(std::get<3>(sizes(ff.hypercubic())) == 12);
	if(ff.set_comm().size() == 2) CHECK(std::get<3>(sizes(ff.hypercubic())) == 6);

	ff.fill(12.2244);

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(ff.matrix()[ii][jj] == 12.2244_a);
		}
	}

	auto zff = complex_field(ff);
	
	static_assert(std::is_same<decltype(zff), basis::field_set<basis::real_space, complex>>::value, "complex() should return a complex field");
		
	CHECK(std::get<1>(sizes(zff.hypercubic())) == 11);
	CHECK(std::get<2>(sizes(zff.hypercubic())) == 20);

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(real(zff.matrix()[ii][jj]) == 12.2244_a);
			CHECK(imag(zff.matrix()[ii][jj]) == 0.0_a);
		}
	}

	auto dff = real_field(zff);

	static_assert(std::is_same<decltype(dff), basis::field_set<basis::real_space, double>>::value, "real() should return a double field");

	CHECK(std::get<1>(sizes(dff.hypercubic())) == 11);
	CHECK(std::get<2>(sizes(dff.hypercubic())) == 20);

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(dff.matrix()[ii][jj] == 12.2244_a);
		}
	}

	dff.fill(0.0);

	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(dff.matrix()[ii][jj] == 0.0_a);
		}
	}
	
	zff.fill(0.0);
	
	for(int ii = 0; ii < ff.basis().part().local_size(); ii++){
		for(int jj = 0; jj < ff.set_part().local_size(); jj++){
			CHECK(real(zff.matrix()[ii][jj]) == 0.0_a);
			CHECK(imag(zff.matrix()[ii][jj]) == 0.0_a);
		}
	}

	basis::field_set<basis::real_space, double> red(basis::field_set<basis::real_space, double>(ff), parallel::cartesian_communicator<2>{comm, {boost::mpi3::fill, 1}});

	CHECK(red.basis().local_size() == red.matrix().size());
	CHECK(red.local_set_size() == (~red.matrix()).size());

	basis::field_set<basis::real_space, double> rr(rs, 12, cart_comm);
	rr.fill(1.0/set_comm.size());
	rr.all_reduce(set_comm);

	for(int ii = 0; ii < rr.basis().part().local_size(); ii++){
		for(int jj = 0; jj < rr.set_part().local_size(); jj++){
			CHECK(rr.matrix()[ii][jj] == 1.0_a);
		}
	}
	
}
#endif
