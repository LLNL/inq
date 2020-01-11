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

#include <basis/field_set.hpp>
#include <cstdlib>

namespace operations {

  template <class matrix_type>
	class matrix_operator {

	public:

    matrix_operator(matrix_type && matrix):
      matrix_(matrix){

      assert(std::get<0>(sizes(matrix_)) == std::get<1>(sizes(matrix_)));
      
    }
      
    template <class field_set_type>
		field_set_type operator()(field_set_type & phi) const {
      
      assert(std::get<0>(sizes(matrix_)) == phi.basis().size());
      assert(std::get<1>(sizes(matrix_)) == phi.basis().size());

      using boost::multi::blas::gemm;
      using boost::multi::blas::hermitized;
    
      field_set_type mphi = phi;
      mphi.matrix() = gemm(phi.basis().volume_element(), matrix_, phi.matrix());
      //mphi.matrix() = gemm(matrix_, phi.matrix());

      return mphi;      
    }

	private:

    matrix_type matrix_;
    
	};
	
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::matrix_operator", "[operations::matrix_operator]") {

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
    
    operations::matrix_operator<math::array<double, 2>> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, double> aa(bas, nvec);

    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2);
      }
    }
    
    auto bb = mo(aa);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        REQUIRE(bb.matrix()[ip][ivec] == Approx((ip + 2.0)*bas.volume_element()*aa.matrix()[ip][ivec]));
      }
    }

  }

  SECTION("Laplacian matrix real"){
  
    math::array<double, 2> matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        matrix[ip][jp] = 0.0;
        if(ip == jp) matrix[ip][jp] = -1.0/bas.volume_element();
        if(ip == jp + 1 or ip == jp - 1) matrix[ip][jp] = 2.0/bas.volume_element();
      }
    }
    
    operations::matrix_operator<math::array<double, 2>> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, double> aa(bas, nvec);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2);
      }
    }    

    auto bb = mo(aa);
    
    for(int ip = 1; ip < npoint - 1; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        REQUIRE(bb.matrix()[ip][ivec] == Approx(2.0*aa.matrix()[ip - 1][ivec] - 1.0*aa.matrix()[ip][ivec] + 2.0*aa.matrix()[ip + 1][ivec]));
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
    
    operations::matrix_operator<math::array<complex, 2>> mo(std::move(matrix));
    
    basis::field_set<basis::trivial, complex> aa(bas, nvec);

    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        aa.matrix()[ip][ivec] = sin((ip + 1.0)*0.5)*cos((ivec + 0.3)*1.2)*exp(complex(0.0, (ip + ivec)*10.0));
      }
    }
    
    auto bb = mo(aa);
    
    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        REQUIRE(real(bb.matrix()[ip][ivec]) == Approx(real(complex(ip + 2.0, 0.3*ip - 6.7)*bas.volume_element()*aa.matrix()[ip][ivec])));
        REQUIRE(imag(bb.matrix()[ip][ivec]) == Approx(imag(complex(ip + 2.0, 0.3*ip - 6.7)*bas.volume_element()*aa.matrix()[ip][ivec])));
      }
    }

  }

  SECTION("Periodic Laplacian matrix complex"){
  
    math::array<complex, 2> matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        matrix[ip][jp] = 0.0;
        if(ip == jp) matrix[ip][jp] = -1.0/bas.volume_element();
        if(ip == jp + 1 or ip == jp - 1) matrix[ip][jp] = 2.0/bas.volume_element();
      }
    }
    //the periodic part
    matrix[0][npoint - 1] = 2.0/bas.volume_element();
    matrix[npoint - 1][0] = 2.0/bas.volume_element();
    
    operations::matrix_operator<math::array<complex, 2>> mo(std::move(matrix));
    
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
        REQUIRE(real(bb.matrix()[ip][ivec]) == Approx(real(2.0*aa.matrix()[prev][ivec] - 1.0*aa.matrix()[ip][ivec] + 2.0*aa.matrix()[next][ivec])));
        REQUIRE(imag(bb.matrix()[ip][ivec]) == Approx(imag(2.0*aa.matrix()[prev][ivec] - 1.0*aa.matrix()[ip][ivec] + 2.0*aa.matrix()[next][ivec])));
      }
    }

  }
  
}


#endif

#endif
