/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__ORTHOGONALIZE
#define INQ__OPERATIONS__ORTHOGONALIZE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>
#include <basis/field_set.hpp>
#include <matrix/cholesky.hpp>
#include <math/complex.hpp>
#include <operations/overlap.hpp>
#include <operations/rotate.hpp>

#include <cstdlib>
#include <utils/profiling.hpp>

namespace inq {
namespace operations {

template <class field_set_type>
void orthogonalize(field_set_type & phi, bool nocheck = false){
	CALI_CXX_MARK_FUNCTION;

	auto olap = overlap(phi);
	matrix::cholesky(olap);
	operations::rotate_trs(olap, phi);
}

}
}
#endif

#ifdef INQ_OPERATIONS_ORTHOGONALIZE_UNIT_TEST
#undef INQ_OPERATIONS_ORTHOGONALIZE_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <operations/randomize.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	systems::box box = systems::box::cubic(6.3_b);
	basis::real_space pw(box, /*spacing =*/ 0.44428829, comm);

	SECTION("Dimension 3"){
		basis::field_set<basis::real_space, complex> phi(pw, 3);
		
		operations::randomize(phi);
		
		operations::orthogonalize(phi);
		
		auto olap = operations::overlap(phi);
		auto olap_array = matrix::all_gather(olap);
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					CHECK(real(olap_array[ii][ii]) == 1.0_a);
					CHECK(fabs(imag(olap_array[ii][ii])) < 1e-14);
			} else {
					CHECK(fabs(olap_array[ii][jj]) < 1e-14);
				}
			}
		}
	}
	
	SECTION("Dimension 2, warning test"){
		basis::field_set<basis::real_space, complex> phi(pw, 2);

		operations::randomize(phi);
		
		int nbasis = std::get<0>(sizes(phi.matrix()));
		
		for(int ii = 0 ; ii < phi.set_size(); ii++){
		  for(int jj = 0; jj < nbasis; jj++){
		    phi.matrix().transposed()[ii][jj]=0.0;
		    //if(ii == jj) phi.matrix().transposed()[ii][jj] = 1.0;
		  }
		}
		phi.matrix().transposed()[0][0] = 1.0;
		phi.matrix().transposed()[1][0] = 1.0;
		
		operations::orthogonalize(phi, /* ncheck = */ true);
		
		//auto olap = operations::overlap(phi);
		
	}

	SECTION("Dimension 100"){
		basis::field_set<basis::real_space, complex> phi(pw, 100);
		
		operations::randomize(phi);
		
		operations::orthogonalize(phi);
		
		auto olap = operations::overlap(phi);
		auto olap_array = matrix::all_gather(olap);
				
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					CHECK(real(olap_array[ii][ii]) == 1.0_a);
					CHECK(fabs(imag(olap_array[ii][ii])) < 1e-14);
				} else {
					CHECK(fabs(olap_array[ii][jj]) < 1e-13);
				}
			}
		}
	}


	SECTION("Dimension 37 - double orthogonalize"){
		basis::field_set<basis::real_space, complex> phi(pw, 37);
		
		operations::randomize(phi);
		
		operations::orthogonalize(phi);
		operations::orthogonalize(phi);
		
		auto olap = operations::overlap(phi);
		auto olap_array = matrix::all_gather(olap);
		
		for(int ii = 0; ii < phi.set_size(); ii++){
			for(int jj = 0; jj < phi.set_size(); jj++){
				if(ii == jj) {
					CHECK(real(olap_array[ii][ii]) == 1.0_a);
					CHECK(fabs(imag(olap_array[ii][ii])) < 1e-16);
				} else {
					CHECK(fabs(olap_array[ii][jj]) < 5e-16);
				}
			}
		}
	}
	
}
#endif
