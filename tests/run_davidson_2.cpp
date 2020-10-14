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



//#include <iostream>
//#include <math/complex.hpp>
#include <eigensolvers/davidson.hpp>
#include <basis/trivial.hpp>
#include <operations/matrix_operator.hpp>
#include <operations/diag_prec.hpp>
#include <multi/adaptors/blas.hpp>
#include <operations/qrfactorize.hpp>
namespace inq{
 
}
int main(int argc, char ** argv){

	using namespace inq;
	
	boost::mpi3::environment env(argc, argv);
	boost::mpi3::communicator comm_world = boost::mpi3::environment::get_world_instance();
			
	const int npoint = 10000;
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
  
	math::array<complex, 2> ham_matrix({npoint, npoint});
    
	for(int ip = 0; ip < npoint; ip++){
	  for(int jp = ip; jp < npoint; jp++){
	    ham_matrix[ip][jp] = (ip + 1)*exp(-0.5*fabs(jp-ip));
	    ham_matrix[jp][ip] = ham_matrix[ip][jp];
	  }
	}
    
	operations::matrix_operator<complex> ham_op(std::move(ham_matrix));
    
	basis::field_set<basis::trivial, complex> phi(bas, nvec, comm_world);
	for(int ip = 0; ip < npoint; ip++){
	  for(int ivec = 0; ivec < nvec; ivec++){
	    phi.matrix()[ip][ivec] = exp(-1.0*abs(ivec-ip));
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

	eigensolvers::davidson(ham_op, identity, phi);

	
}
