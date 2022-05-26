/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__ROTATE
#define OPERATIONS__ROTATE

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa.

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

#include <inq_config.h>

#include <math/array.hpp>
#include <utils/profiling.hpp>

#include <cassert>

namespace inq {
namespace operations {

//////////////////////////////////////////////////////////////////////////////////

template <class MatrixType, class FieldSetType>
void rotate(MatrixType const & rotation, FieldSetType & phi){
	
	CALI_CXX_MARK_SCOPE("operations::rotate");

	namespace blas = boost::multi::blas;
	//OPTIMIZATION: here we don't need to make a full copy.
	phi.matrix() = blas::gemm(1., phi.matrix(), blas::H(rotation.array()));
}

//////////////////////////////////////////////////////////////////////////////////

template <class MatrixType, class FieldSetType>
void rotate_trs(MatrixType const & rotation, FieldSetType & phi){
	
	CALI_CXX_MARK_SCOPE("operations::rotate_trs");

	namespace blas = boost::multi::blas;
	blas::trsm(blas::side::right, blas::filling::upper, 1.0, blas::H(rotation.array()), phi.matrix());
}

//////////////////////////////////////////////////////////////////////////////////

}
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_ROTATE_UNIT_TEST
#undef INQ_OPERATIONS_ROTATE_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>
#include <math/subspace_matrix.hpp>
#include <operations/orthogonalize.hpp>
#include <operations/randomize.hpp>

TEST_CASE("function operations::rotate", "[operations::rotate]") {
	
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 100;
	const int nvec = 4;
			
	auto comm = boost::mpi3::environment::get_world_instance();

	auto parstates = comm.size();
	if(comm.size() >= 1) parstates = 1;
	
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {parstates, boost::mpi3::fill});
	auto basis_comm = cart_comm.axis(1);
	
	basis::trivial bas(npoint, basis_comm);
	
	SECTION("rotate double"){
		
		math::subspace_matrix<double> rot(cart_comm, nvec, 0.0);
		
		for(int ii = 0; ii < nvec; ii++){
			for(int jj = 0; jj < nvec; jj++){
				rot.array()[ii][jj] = ii + 1;
			}
		}

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
	
	SECTION("rotate complex"){
		
		math::subspace_matrix<complex> rot(cart_comm, nvec, 0.0);
		
		for(int ii = 0; ii < nvec; ii++){
			for(int jj = 0; jj < nvec; jj++){
				rot.array()[ii][jj] = ii + 1;
			}
		}
		
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
	
	
}

#endif
#endif
