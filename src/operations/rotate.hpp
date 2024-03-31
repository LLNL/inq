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
			auto res = ncclReduce(raw_pointer_cast(block.data_elements()), raw_pointer_cast(phi_matrix.base()),
														block.num_elements()*sizeof(typename PhiMatrix::element_type)/sizeof(double), ncclDouble, ncclSum, istep, &set_comm.nccl_comm(), 0);		 
			assert(res == ncclSuccess);
			gpu::sync();
#else
			set_comm.reduce_n(raw_pointer_cast(block.data_elements()), block.num_elements(), raw_pointer_cast(phi_matrix.base()), std::plus<>{},  /* root = */ istep);
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

template <class MatrixType, class SetComm, class PhiSetPart, class PhiMatrix, class RotPhiMatrix>
void rotate_impl(MatrixType const & rotation, SetComm & set_comm,
								 PhiSetPart const phi_set_part, PhiMatrix const & phi_matrix, RotPhiMatrix & rotphi_matrix,
								 typename PhiMatrix::element_type const & alpha, typename PhiMatrix::element_type const & beta){

	namespace blas = boost::multi::blas;

	CALI_CXX_MARK_SCOPE("operations::rotate(5arg)");

	auto rotation_array = matrix::all_gather(rotation);
	
	if(not phi_set_part.parallel()){

		blas::gemm(alpha, phi_matrix, rotation_array, beta, rotphi_matrix);
		
	} else {

		set_comm.nccl_init();
		
		for(int istep = 0; istep < set_comm.size(); istep++){
				
			auto block = +blas::gemm(alpha, phi_matrix, rotation_array({phi_set_part.start(), phi_set_part.end()}, {phi_set_part.start(istep), phi_set_part.end(istep)}));
			
			assert(block.extensions() == phi_matrix.extensions());

			if(istep == set_comm.rank() and beta != 0.0){
				gpu::run(phi_set_part.local_size(), phi_matrix.size(),
								 [blo = begin(block), rot = begin(rotphi_matrix), beta] GPU_LAMBDA (auto ist, auto ip){
									 blo[ip][ist] += beta*rot[ip][ist];
								 });
			}
			
			CALI_CXX_MARK_SCOPE("operations::rotate(5arg)_reduce");
#ifdef ENABLE_NCCL
			auto res = ncclReduce(raw_pointer_cast(block.data_elements()), raw_pointer_cast(rotphi_matrix.base()),
														block.num_elements()*sizeof(typename RotPhiMatrix::element_type)/sizeof(double), ncclDouble, ncclSum, istep, &set_comm.nccl_comm(), 0);		 
			assert(res == ncclSuccess);
			gpu::sync();			
#else
			set_comm.reduce_n(raw_pointer_cast(block.data_elements()), block.num_elements(), raw_pointer_cast(rotphi_matrix.base()), std::plus<>{}, /* root = */ istep);
#endif
		}
	}
	
}

//////////////////////////////////////////////////////////////////////////////////

template <class Matrix, class Basis, class Type, class ScalType>
void rotate(Matrix const & rotation, basis::field_set<Basis, Type> const & phi, basis::field_set<Basis, Type> & rotphi, ScalType const & alpha = 1.0, ScalType const & beta = 0.0){
	rotate_impl(rotation, phi.set_comm(), phi.set_part(), phi.matrix(), rotphi.matrix(), alpha, beta);
}

//////////////////////////////////////////////////////////////////////////////////

template <class Matrix, class Basis, class Type, class ScalType>
void rotate(Matrix const & rotation, states::orbital_set<Basis, Type> const & phi, states::orbital_set<Basis, Type> & rotphi, ScalType const & alpha = 1.0, ScalType const & beta = 0.0){
	auto rotphi_matrix =  rotphi.basis_spinor_matrix();
	rotate_impl(rotation, phi.set_comm(), phi.spinor_set_part(), phi.basis_spinor_matrix(), rotphi_matrix, alpha, beta);
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
		
		SECTION("rotate orbital_set spinor double"){
			
			gpu::array<double, 2> rot_array({nvec, nvec}, 0.0);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = ii + 1;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, double> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto jjg = aa.spinor_set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.spinor_matrix()[ip][0][jj] = (jjg.value() + 1.0)*(ipg.value() + 1);
					aa.spinor_matrix()[ip][1][jj] = -2.0*(jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto jjg = aa.spinor_set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					CHECK(aa.spinor_matrix()[ip][0][jj] == (ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
					CHECK(aa.spinor_matrix()[ip][1][jj] == -2.0*(ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
				}
			}
		
		}

		SECTION("rotate orbital_set spinor complex"){
		
			gpu::array<complex, 2> rot_array({nvec, nvec}, 0.0);
		
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = ii + 1;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto jjg = aa.spinor_set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.spinor_matrix()[ip][0][jj] = complex{0.0, (jjg.value() + 1.0)*(ipg.value() + 1)};
					aa.spinor_matrix()[ip][1][jj] = complex{-2.0*(jjg.value() + 1.0)*(ipg.value() + 1), 0.0};
				}
			}

			operations::rotate(rot, aa);

			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto jjg = aa.spinor_set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					CHECK(real(aa.spinor_matrix()[ip][0][jj]) == 0.0);
					CHECK(imag(aa.spinor_matrix()[ip][0][jj]) == (ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
					CHECK(real(aa.spinor_matrix()[ip][1][jj]) == -2.0*(ipg.value() + 1.0)*(jjg.value() + 1.0)*nvec*(nvec + 1.0)/2.0);
					CHECK(imag(aa.spinor_matrix()[ip][1][jj]) == 0.0);
				}
			}
		
		}

		SECTION("rotate_trs orbital_set spinor double"){
			
			gpu::array<double, 2> rot_array({nvec, nvec}, 0.0);;
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = 1.0;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, double> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto jjg = aa.spinor_set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.spinor_matrix()[ip][0][jj] = (jjg.value() + 1.0)*(ipg.value() + 1);
					aa.spinor_matrix()[ip][1][jj] = -0.5*(jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate_trs(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto ipg = bas.part().local_to_global(ip);
					CHECK(aa.spinor_matrix()[ip][0][jj] == (ipg.value() + 1.0));
					CHECK(aa.spinor_matrix()[ip][1][jj] == -0.5*(ipg.value() + 1.0));
				}
			}
			
		}

		SECTION("rotate_trs orbital_set spinor complex"){
			
			gpu::array<complex, 2> rot_array({nvec, nvec}, 0.0);
			
			for(int ii = 0; ii < nvec; ii++){
				for(int jj = 0; jj < nvec; jj++){
					rot_array[ii][jj] = 1.0;
				}
			}

			auto rot = matrix::scatter(cart_comm, rot_array, /* root = */ 0);
			
			states::orbital_set<basis::trivial, complex> aa(bas, nvec, /*spinor_dim = */ 2, /*kpoint = */ vector3<double, covariant>{0.0, 0.0, 0.0}, /*spin_index = */ 0, cart_comm);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto jjg = aa.spinor_set_part().local_to_global(jj);
					auto ipg = bas.part().local_to_global(ip);
					aa.spinor_matrix()[ip][0][jj] = complex{0.0, (jjg.value() + 1.0)*(ipg.value() + 1)};
					aa.spinor_matrix()[ip][1][jj] = 0.4*(jjg.value() + 1.0)*(ipg.value() + 1);
				}
			}

			operations::rotate_trs(rot, aa);
			
			for(int ip = 0; ip < bas.part().local_size(); ip++){
				for(int jj = 0; jj < aa.local_spinor_set_size(); jj++){
					auto ipg = bas.part().local_to_global(ip);
					CHECK(real(aa.spinor_matrix()[ip][0][jj]) == Approx(0.0));
					CHECK(imag(aa.spinor_matrix()[ip][0][jj]) == Approx((ipg.value() + 1.0)));
					CHECK(real(aa.spinor_matrix()[ip][1][jj]) == Approx(0.4*(ipg.value() + 1.0)));
					CHECK(imag(aa.spinor_matrix()[ip][1][jj]) == Approx(0.0));
				}
			}
			
		}


	}


}
#endif
