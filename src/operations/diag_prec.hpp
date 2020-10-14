/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__DIAG_PREC
#define INQ__OPERATIONS__DIAG_PREC

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <operations/space.hpp>

#include <cstdlib>

namespace inq {
namespace operations {

  class diag_prec {
	
    //Implements a simple diagonal preconditioner
  public:
			
    template <class vectype>
    diag_prec(vectype && diagH) {
      int ndim = std::get<0>(sizes(diagH));
    }
    template <class vectype, class bas_type, class data_type>
    void apply_prec(vectype && eigW, basis::field_set<bas_type, data_type> & phi){
      int nvec = std::get<0>(sizes(eigW));	
      int ndim = std::get<1>(sizes(phi.matrix()));
    }

private:

};
	
}
}

#ifdef INQ_OPERATIONS_DIAG_PREC_UNIT_TEST
#undef INQ_OPERATIONS_DIAG_PREC_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("function operations::diag_prec", "[diag_prec]") {

	using namespace inq;
	using namespace Catch::literals;

}


#endif

#endif
