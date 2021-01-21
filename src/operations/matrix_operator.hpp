/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__MATRIX_OPERATOR
#define OPERATIONS__MATRIX_OPERATOR

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

#include <math/array.hpp>

#include <cstdlib>

namespace inq {
namespace operations {

template <class type>
class matrix_operator {
	
public:

	matrix_operator(math::array<type, 2> && matrix):
		matrix_(matrix){

		assert(std::get<0>(sizes(matrix_)) == std::get<1>(sizes(matrix_)));
      
	}
      
	template <class field_set_type>
	field_set_type operator()(const field_set_type & phi) const {
      
		assert(std::get<0>(sizes(matrix_)) == phi.basis().size());
		assert(std::get<1>(sizes(matrix_)) == phi.basis().size());

		using boost::multi::blas::gemm;
		using boost::multi::blas::hermitized;
    
		field_set_type mphi = phi;
		gemm(1.0, matrix_, phi.matrix(),	0.0, mphi.matrix());

		return mphi;      
	}

private:

	math::array<type, 2> matrix_;
    
};

}
}

#ifdef INQ_OPERATIONS_MATRIX_OPERATOR_UNIT_TEST
#undef INQ_OPERATIONS_MATRIX_OPERATOR_UNIT_TEST

#include <basis/trivial.hpp>
#include <basis/field_set.hpp>
#include <catch2/catch.hpp>

TEST_CASE("function operations::matrix_operator", "[operations::matrix_operator]") {

	using namespace inq;
	using namespace Catch::literals;

  const int npoint = 100;
  const int nvec = 12;
  
  basis::trivial bas(npoint);

  SECTION("Diagonal matrix real"){

    math::array<double, 2> matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        matrix[ip][jp] = 0.0;
        if(ip == jp) matrix[ip][jp] = ip + 2.0;
      }
    }
    
    operations::matrix_operator<double> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, double> aa(bas, nvec);

    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2);
      }
    }
    
    auto bb = mo(aa);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        CHECK(bb.matrix()[ip][ivec] == Approx((ip + 2.0)*aa.matrix()[ip][ivec]));
      }
    }

  }

  SECTION("Laplacian matrix real"){
  
    math::array<double, 2> matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        matrix[ip][jp] = 0.0;
        if(ip == jp) matrix[ip][jp] = -1.0;
        if(ip == jp + 1 or ip == jp - 1) matrix[ip][jp] = 2.0;
      }
    }
    
    operations::matrix_operator<double> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, double> aa(bas, nvec);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2);
      }
    }    

    auto bb = mo(aa);
    
    for(int ip = 1; ip < npoint - 1; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        CHECK(bb.matrix()[ip][ivec] == Approx(2.0*aa.matrix()[ip - 1][ivec] - 1.0*aa.matrix()[ip][ivec] + 2.0*aa.matrix()[ip + 1][ivec]));
      }
    }

  }

  SECTION("Diagonal matrix complex"){

    math::array<complex, 2> matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        matrix[ip][jp] = 0.0;
        if(ip == jp) matrix[ip][jp] = complex(ip + 2.0, 0.3*ip - 6.7);
      }
    }
    
    operations::matrix_operator<complex> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, complex> aa(bas, nvec);

    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2)*exp(complex(0.0, (ip + ivec)*10.0));
      }
    }
    
    auto bb = mo(aa);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        CHECK(real(bb.matrix()[ip][ivec]) == Approx(real(complex(ip + 2.0, 0.3*ip - 6.7)*aa.matrix()[ip][ivec])));
        CHECK(imag(bb.matrix()[ip][ivec]) == Approx(imag(complex(ip + 2.0, 0.3*ip - 6.7)*aa.matrix()[ip][ivec])));
      }
    }

  }

  SECTION("Periodic Laplacian matrix complex"){
  
    math::array<complex, 2> matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        matrix[ip][jp] = 0.0;
        if(ip == jp) matrix[ip][jp] = -1.0;
        if(ip == jp + 1 or ip == jp - 1) matrix[ip][jp] = 2.0;
      }
    }
    //the periodic part
    matrix[0][npoint - 1] = 2.0;
    matrix[npoint - 1][0] = 2.0;
    
    operations::matrix_operator<complex> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, complex> aa(bas, nvec);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2)*exp(complex(0.0, (ip + ivec)*10.0));
      }
    }    

    auto bb = mo(aa);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        auto prev = ip - 1;
        auto next = ip + 1;
        if(prev == -1) prev = npoint - 1;
        if(next == npoint) next = 0;
        CHECK(real(bb.matrix()[ip][ivec]) == Approx(real(2.0*aa.matrix()[prev][ivec] - 1.0*aa.matrix()[ip][ivec] + 2.0*aa.matrix()[next][ivec])));
        CHECK(imag(bb.matrix()[ip][ivec]) == Approx(imag(2.0*aa.matrix()[prev][ivec] - 1.0*aa.matrix()[ip][ivec] + 2.0*aa.matrix()[next][ivec])));
      }
    }

  }
  
}


#endif

#endif
