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

template <class MatrixType, class FieldSetType>
void rotate(MatrixType const & rotation, FieldSetType & phi){
	
	CALI_CXX_MARK_SCOPE("operations::rotate");

	namespace blas = boost::multi::blas;

	//OPTIMIZATION: here we don't need to make a full copy.
	phi.matrix() = blas::gemm(1., phi.matrix(), blas::H(rotation.array()));
	
}

}
}

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

}

#endif
#endif
