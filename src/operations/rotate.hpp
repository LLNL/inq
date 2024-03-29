/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__ROTATE
#define OPERATIONS__ROTATE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/copy.hpp>
#include <gpu/array.hpp>
#include <matrix/invert_triangular.hpp>
#include <utils/profiling.hpp>

#include <cassert>

namespace inq {
namespace operations {

//////////////////////////////////////////////////////////////////////////////////

template <class MatrixType, class SetPart, class SetComm, class PhiMatrix>
void rotate_impl(MatrixType const & rotation, SetPart const & set_part, SetComm & set_comm, PhiMatrix & phi_matrix){
	
	CALI_CXX_MARK_SCOPE("operations::rotate(2arg)");

	namespace blas = boost::multi::blas;

	auto rotation_array = matrix::all_gather(rotation);
	
	if(set_part.parallel()){

		set_comm.nccl_init();

		// The direct copy is slow with multi right now: auto copy = phi_matrix;
		gpu::array<typename PhiMatrix::element_type, 2> copy({phi_matrix.size(), set_part.local_size()});
		gpu::copy(phi_matrix.size(), set_part.local_size(), phi_matrix, copy);
		
		for(int istep = 0; istep < set_part.comm_size(); istep++){
			
			auto block = +blas::gemm(1.0, copy, blas::H(rotation_array({set_part.start(istep), set_part.end(istep)}, {set_part.start(), set_part.end()}))); 
			
			assert(block.extensions() == phi_matrix.extensions());

			CALI_CXX_MARK_SCOPE("operations::rotate(2arg)_reduce");
#ifdef ENABLE_NCCL
			auto res = ncclReduce(raw_pointer_cast(block.data_elements()), raw_pointer_cast(phi_matrix.data_elements()),
														block.num_elements()*sizeof(typename FieldSetType::element_type)/sizeof(double), ncclDouble, ncclSum, istep, &set_comm.nccl_comm(), 0);		 
			assert(res == ncclSuccess);
			gpu::sync();
#else
			set_comm.reduce_n(raw_pointer_cast(block.data_elements()), block.num_elements(), raw_pointer_cast(phi_matrix.base()), std::plus<>{}, istep);
#endif
		}

	} else {
		phi_matrix = +blas::gemm(1.0, phi_matrix, blas::H(rotation_array));
	}
	
}

//////////////////////////////////////////////////////////////////////////////////

template <class Matrix, class Basis, class Type>
void rotate(Matrix const & rotation, basis::field_set<Basis, Type> & phi){
	rotate_impl(rotation, phi.set_part(), phi.set_comm(), phi.matrix());
}

//////////////////////////////////////////////////////////////////////////////////

template <class Matrix, class Basis, class Type>
void rotate(Matrix const & rotation, states::orbital_set<Basis, Type> & phi){
	auto phimatrix = phi.basis_spinor_matrix();
	rotate_impl(rotation, phi.spinor_set_part(), phi.set_comm(), phimatrix);
}

//////////////////////////////////////////////////////////////////////////////////

template <class MatrixType, class FieldSetType1, class FieldSetType2>
void rotate(MatrixType const & rotation, FieldSetType1 const & phi, FieldSetType2 & rotphi, typename FieldSetType1::element_type const & alpha = 1.0, typename FieldSetType1::element_type const & beta = 0.0){
	namespace blas = boost::multi::blas;

	CALI_CXX_MARK_SCOPE("operations::rotate(5arg)");

	auto rotation_array = matrix::all_gather(rotation);
	
	if(not phi.set_part().parallel()){

		blas::gemm(alpha, phi.matrix(), rotation_array, beta, rotphi.matrix());
		
	} else {

		phi.set_comm().nccl_init();
		
		for(int istep = 0; istep < phi.set_part().comm_size(); istep++){
				
			auto block = +blas::gemm(alpha, phi.matrix(), rotation_array({phi.set_part().start(), phi.set_part().end()}, {phi.set_part().start(istep), phi.set_part().end(istep)}));
			
			assert(block.extensions() == phi.matrix().extensions());

			if(istep == rotphi.set_comm().rank() and beta != 0.0){
				gpu::run(phi.local_set_size(), phi.basis().local_size(),
								 [blo = begin(block), rot = begin(rotphi.matrix()), beta] GPU_LAMBDA (auto ist, auto ip){
									 blo[ip][ist] += beta*rot[ip][ist];
								 });
			}
			
			CALI_CXX_MARK_SCOPE("operations::rotate(5arg)_reduce");
#ifdef ENABLE_NCCL
			auto res = ncclReduce(raw_pointer_cast(block.data_elements()), raw_pointer_cast(rotphi.matrix().data_elements()),
														block.num_elements()*sizeof(typename FieldSetType1::element_type)/sizeof(double), ncclDouble, ncclSum, istep, &phi.set_comm().nccl_comm(), 0);		 
			assert(res == ncclSuccess);
			gpu::sync();			
#else
			phi.set_comm().reduce_n(raw_pointer_cast(block.data_elements()), block.num_elements(), raw_pointer_cast(rotphi.matrix().data_elements()), std::plus<>{}, istep);
#endif
		}
	}
	
}

//////////////////////////////////////////////////////////////////////////////////

template <class MatrixType, class FieldSetType>
void rotate_trs(MatrixType const & rotation, FieldSetType & phi){
	
	CALI_CXX_MARK_SCOPE("operations::rotate_trs");

	auto invrot = rotation;
	
	namespace blas = boost::multi::blas;
	matrix::invert_triangular(invrot);
	rotate(invrot, phi);
}

//////////////////////////////////////////////////////////////////////////////////

}
}
#endif

#ifdef INQ_OPERATIONS_ROTATE_UNIT_TEST
#undef INQ_OPERATIONS_ROTATE_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>
#include <operations/orthogonalize.hpp>
#include <operations/randomize.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 100;
	const int nvec = 4;
			
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	{
		auto parstates = comm.size();
		if(comm.size() == 3 or comm.size() >= 5) parstates = 1;
		if(comm.size() == 2 or comm.size() == 4) parstates = 2;
		
		parallel::cartesian_communicator<2> cart_comm(comm, {boost::mpi3::fill, parstates});
		auto basis_comm = basis::basis_subcomm(cart_comm);
		
		basis::trivial bas(npoint, basis_comm);

		SECTION("rotate field_set double"){
			
			gpu::array<double, 2> rot_array({nvec, nvec}, 0.0);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = ii + 1;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = (jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					CHECK(aa.matrix()[ip][jj] == (ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
				}
			}
		
		}
		
		SECTION("rotate field_set complex"){
		
			gpu::array<complex, 2> rot_array({nvec, nvec}, 0.0);
		
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = ii + 1;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = complex{0.0, (jjg.value() + 1.0)*(ipg.value() + 1)};
				}
			}

			operations::rotate(rot, aa);

			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					CHECK(real(aa.matrix()[ip][jj]) == 0.0);
					CHECK(imag(aa.matrix()[ip][jj]) == (ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
				}
			}
		
		}
		
		SECTION("rotate_trs field_set double"){
			
			gpu::array<double, 2> rot_array({nvec, nvec}, 0.0);;
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = 1.0;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = (jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate_trs(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto ipg = bas.part().local_to_global(ip);
					CHECK(aa.matrix()[ip][jj] == (ipg.value() + 1.0));
				}
			}
			
		}
		
		SECTION("rotate_trs field_set complex"){
			
			gpu::array<complex, 2> rot_array({nvec, nvec}, 0.0);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = 1.0;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = complex{0.0, (jjg.value() + 1.0)*(ipg.value() + 1)};
				}
			}

			operations::rotate_trs(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto ipg = bas.part().local_to_global(ip);
					CHECK(real(aa.matrix()[ip][jj]) == 0.0);
					CHECK(imag(aa.matrix()[ip][jj]) == (ipg.value() + 1.0));
				}
			}
			
		}

		SECTION("rotate orbital_set double"){
			
			gpu::array<double, 2> rot_array({nvec, nvec}, 0.0);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = ii + 1;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, double> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = (jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					CHECK(aa.matrix()[ip][jj] == (ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
				}
			}
		
		}
		
		SECTION("rotate orbital_set complex"){
		
			gpu::array<complex, 2> rot_array({nvec, nvec}, 0.0);
		
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = ii + 1;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = complex{0.0, (jjg.value() + 1.0)*(ipg.value() + 1)};
				}
			}

			operations::rotate(rot, aa);

			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					CHECK(real(aa.matrix()[ip][jj]) == 0.0);
					CHECK(imag(aa.matrix()[ip][jj]) == (ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
				}
			}
		
		}
		
		SECTION("rotate_trs orbital_set double"){
			
			gpu::array<double, 2> rot_array({nvec, nvec}, 0.0);;
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = 1.0;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, double> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = (jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate_trs(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto ipg = bas.part().local_to_global(ip);
					CHECK(aa.matrix()[ip][jj] == (ipg.value() + 1.0));
				}
			}
			
		}
		
		SECTION("rotate_trs orbital_set complex"){
			
			gpu::array<complex, 2> rot_array({nvec, nvec}, 0.0);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = 1.0;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 1, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto jjg = aa.set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.matrix()[ip][jj] = complex{0.0, (jjg.value() + 1.0)*(ipg.value() + 1)};
				}
			}

			operations::rotate_trs(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_set_size(); jj++){
					auto ipg = bas.part().local_to_global(ip);
					CHECK(real(aa.matrix()[ip][jj]) == 0.0);
					CHECK(imag(aa.matrix()[ip][jj]) == (ipg.value() + 1.0));
				}
			}
			
		}
		
	}
	
}
#endif
