/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019-2020 Xavier Andrade

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


//#include <systems/ions.hpp>
//#include <systems/electrons.hpp>
//#include <config/path.hpp>
//#include <input/atom.hpp>
//#include <utils/match.hpp>
//#include <iostream>
//#include <math/complex.hpp>
//#include <multi/adaptors/blas.hpp>
#include <eigensolvers/davidson.hpp>
#include <basis/trivial.hpp>
#include <operations/matrix_operator.hpp>
#include <operations/diag_prec.hpp>
//#include <multi/adaptors/totalview.hpp> 
#include <multi/adaptors/blas.hpp>
#include <operations/qrfactorize.hpp>
#include <operations/overlap.hpp>
namespace inq{
  // template <class array_type>
  // auto print_matrix(array_type & A){
  //   int nrowA=std::get<0>(sizes(A));
  //   int ncolA=std::get<1>(sizes(A));
  //   for(int i=0; i<nrowA ; i++){
  //     for(int j=0; j<ncolA ; j++){
  // 	std::printf("%4.1f ",A[i][j]);
  //     }
  //     std::cout << '\n';
  //   }      
  // }
  // template <class array_type>
  // auto print_vec(array_type & A){
  //   int nrowA=std::get<0>(sizes(A));
  //   for(int i=0; i<nrowA ; i++){
  //     std::printf("%10.6f ",A[i]);
  //   }
  //   std::cout << '\n';
  // }
}
int main(int argc, char ** argv){

	using namespace inq;
	
	inq::input::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
			
	const int npoint = 1000;
	const int nvec = 9;
  
	basis::trivial bas(npoint);
	
	math::array<complex, 2> identity_matrix({npoint, npoint});
  
	for(int ip = 0; ip < npoint; ip++){
	  for(int jp = 0; jp < npoint; jp++){
	    identity_matrix[ip][jp] = 0.0;
	    if(ip == jp) identity_matrix[ip][jp] = 1.0;
	  }
	}
  
	operations::matrix_operator<complex> identity(std::move(identity_matrix));
//   SECTION("Diagonal matrix complex"){
  
	// math::array<complex, 2> diagonal_matrix({npoint, npoint});
    
	// for(int ip = 0; ip < npoint; ip++){
	//   for(int jp = 0; jp < npoint; jp++){
	//     diagonal_matrix[ip][jp] = 0.0;
	//     if(ip == jp) diagonal_matrix[ip][jp] = ip + 1.0;
	//   }
	// }
    
	// operations::matrix_operator<complex> diagonal_op(std::move(diagonal_matrix));
	//  SECTION("Diagonal matrix complex"){
  
	math::array<complex, 2> ham_matrix({npoint, npoint});
    
	for(int ip = 0; ip < npoint; ip++){
	  for(int jp = ip; jp < npoint; jp++){
	    ham_matrix[ip][jp] = (ip + 1)*exp(-1.0*fabs(jp-ip));
	    ham_matrix[jp][ip] = ham_matrix[ip][jp];
	  }
	}
    
	operations::matrix_operator<complex> ham_op(std::move(ham_matrix));
    
	basis::field_set<basis::trivial, complex> phi(bas, nvec, comm_world);

	for(int ip = 0; ip < npoint; ip++){
	  for(int ivec = 0; ivec < nvec; ivec++){
	    phi.matrix()[ip][ivec] = 0.0;
	    if(ivec == ip){
	      phi.matrix()[ip][ivec] = 1.0;
	    }
	  }
	}

	operations::qrfactorize(phi);
	//normalize phi
	for(int ivec=0; ivec < nvec; ivec++){
	  auto vec = phi.matrix().rotated()[ivec];
	  double nrm;
      	  auto x = boost::multi::blas::nrm2(vec,nrm);
	  phi.matrix().rotated()[ivec] = boost::multi::blas::scal(1.0e0/nrm,vec); ;
	}
	//print_matrix(phi.matrix());

	// basis::field_set<decltype(phi.basis()), complex> phi2(phi.basis(), nvec, phi.full_comm());

	//basis::field_set<decltype(phi.basis()), complex> phi2(bas, nvec);
	//phi2.matrix() = phi.matrix();
	//operations::orthogonalize(phi2);
	//print_matrix(phi2.matrix());
//     const int num_iter = 100;
		
//     for(int iter = 0; iter < num_iter; iter++){
	// math::array<complex, 2> Am = {{1., 0.36787944},{0.36787944, 2.0}};
	
	// auto ae = operations::diagonalize(Am);
	// printf("%16.8f , %16.8f\n", ae[0],ae[1]);

	eigensolvers::davidson(ham_op, identity, phi);

	
}
