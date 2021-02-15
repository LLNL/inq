/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__EIGENSOLVERS__DAVIDSON
#define INQ__EIGENSOLVERS__DAVIDSON

/*
  Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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

#include <math/complex.hpp>
#include <math/vector3.hpp>
#include <math/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/shift.hpp>
//#include <operations/qrfactorize.hpp>
#include <operations/orthogonalize.hpp>
#include <operations/diagonalize.hpp>
#include <operations/overlap.hpp>
#include <operations/overlap_diagonal.hpp>

template <class array_type>
auto print_matrix(array_type & A){
  int nrowA=std::get<0>(sizes(A));
  int ncolA=std::get<1>(sizes(A));
  for(int i=0; i<nrowA ; i++){
    for(int j=0; j<ncolA ; j++){
      std::printf("%10.6f ",A[i][j]);
    }
    std::cout << '\n';
  }      
}
template <class array_type>
auto print_vec(array_type & A){
  int nrowA=std::get<0>(sizes(A));
  for(int i=0; i<nrowA ; i++){
    std::printf("%10.6f ",A[i]);
  }
  std::cout << '\n';
}

namespace inq {
  namespace eigensolvers {

    template <class operator_type, class preconditioner_type, class field_set_type>
    void davidson(const operator_type & ham, const preconditioner_type & prec, field_set_type & phi){
      //Get some sizes
      const int nbasis  = std::get<0>(sizes(phi.matrix()));
      const int nvec = std::get<1>(sizes(phi.matrix()));

      const int naux = 4 * nvec ; //naux is maximum subspace size
      const int nsize = naux + 2* nvec ; //nsize is maximum size of work wavefunction

      //Initialze some book-keeping variables
      int nbase = nvec ; //current subspace size
      int nevf  = 0 ; //number of eigenvectors found
      int nleft = nvec ; // number of eigenvectors left to converge
      
      //allocate auxiliary wave-functions storage
      field_set_type aux_phi(phi.basis(), nsize, phi.full_comm());
      aux_phi.matrix()({0,nbasis},{0,nvec})=phi.matrix()({0,nbasis},{0,nvec});
      math::array<double, 1> eigW(nvec,0); //eigenvalues saved here
      
      const int num_steps = 20;
      namespace blas = boost::multi::blas ;

      for(int istep = 0; istep < num_steps; istep++){

      	//allocate work wave funcs 
      	field_set_type wphi(phi.basis(), nbase, phi.full_comm());
	wphi.matrix()({0,nbasis},{0,nbase})=aux_phi.matrix()({0,nbasis},{0,nbase});

      	auto Wk = ham(wphi);
      	auto Yk = operations::overlap(wphi, Wk);  //subspace hamiltonian
      	auto lk = operations::diagonalize(Yk);
	//print_vec(lk);
		namespace blas = boost::multi::blas;
      	//Ritz vectors
      	auto Xk = blas::gemm(wphi.matrix(), blas::hermitized(Yk)); //Yk now contains subspace eigenvectors
      	//Rotated hphi
      	auto Tk = blas::gemm(Wk.matrix(), blas::hermitized(Yk));
      	//Residuals
      	nleft = nvec - nevf; //update nleft
      	for(int i=0; i<nleft ; i++){
      	  wphi.matrix().rotated()[i] = blas::axpy(lk[i],Xk.rotated()[i], blas::scal(-1.0,Tk.rotated()[i]));   //wphi=l*psi+(-1.0)H*psi
      	}
	
      	//Deflate here
      	int nc=0;
	int notconv = nleft;

      	if(notconv < 0.25*nvec){
	  //print_vec(eigW);
      	  return ;
      	}
	else if( istep == num_steps -1 ){
     		for(int i=0; i<notconv ; i++){
      	    		phi.matrix().rotated()[nevf+i] = Xk.rotated()[i];
      	    		eigW(nevf+i) = lk[i];
		}
		//print_vec(eigW);
		return ;
	} 
      	//correction vectors
      	field_set_type tphi(phi.basis(), notconv, phi.full_comm());

	tphi.matrix()({0,nbasis},{0,notconv})=wphi.matrix()({0,nbasis},{0,notconv});

      	prec(tphi);
      	if(nevf > 0){
      	  aux_phi.matrix()({0,nbasis},{nevf,nevf+nbase})=aux_phi.matrix()({0,nbasis},{0,nbase});
      	  aux_phi.matrix()({0,nbasis},{0,nevf})=phi.matrix()({0,nbasis},{0,nevf});
      	}
      	aux_phi.matrix()({0,nbasis},{nevf+nbase,nevf+nbase+notconv})=tphi.matrix()({0,nbasis},{0,notconv});
      	//orthogonalize everything   //allocate work wave funcs 
        field_set_type ophi(phi.basis(), nevf+nbase+notconv, phi.full_comm());
#ifdef CUDA_ENABLED
	cudaDeviceSynchronize();
#endif
      	ophi.matrix() = aux_phi.matrix()({0, nbasis},{0, nevf+nbase+notconv});
//assert(size(ophi.matrix()) == nbasis);
//assert(size(~ophi.matrix()) == nevf+nbase+notconv);
//assert( ophi.matrix() == ophi.matrix() );
//assert( ophi.matrix()({0, nbasis}, {0, nevf+nbase+notconv}) == ophi.matrix()({0, nbasis}, {0, nevf+nbase+notconv}) );
//assert( ophi.matrix()({0, nbasis}, {0, nevf+nbase+notconv}) == ophi.matrix() );
	//operations::qrfactorize(ophi);
	operations::orthogonalize(ophi, /* nocheck = */ true);
	
	auto qrfnrm = operations::overlap_diagonal(ophi,ophi);

      	double droptol=1.0e-4;
      	int nkeep = 0;
      	//filter by norm of correction vectors
      	for(int i=nevf;i<nevf+nbase+notconv;i++){
      	  if(abs(qrfnrm[i]) > droptol){
      	    aux_phi.matrix().rotated()[nkeep]=ophi.matrix().rotated()[i];
      	    nkeep = nkeep + 1;
      	  }
      	}
	
      	if(nbase + notconv > naux || nc > 0 || nkeep == nbase){
      	  if(nevf > 0){
      	    aux_phi.matrix()({0,nbasis},{0,nevf}) = phi.matrix()({0,nbasis},{0,nevf});
      	  }
      	  nleft = nvec - nevf ;
      	  for(int i=0;i<nleft;i++){
      	    aux_phi.matrix().rotated()[nevf+i]=Xk.rotated()[i];
      	  }
      	  for(int i=0;i<notconv;i++){
      	    aux_phi.matrix().rotated()[nevf+nleft+i]=tphi.matrix().rotated()[i];
      	  }
      	  field_set_type qphi(phi.basis(), nevf+nleft+notconv, phi.full_comm());
      	  qphi.matrix() = aux_phi.matrix()({0,nbasis},{0,nevf+nleft+notconv});
      	  //operations::qrfactorize(qphi);
	  operations::orthogonalize(qphi, /* nocheck = */ true);
	  auto qrfnrm = operations::overlap_diagonal(qphi,qphi);
      	  double droptol=1.0e-4;
      	  int nkeep = 0;
      	  //filter by norm of correction vectors
      	  for(int i=nevf;i<nevf+nleft+notconv;i++){
      	    if(abs(qrfnrm[i]) > droptol){
      	      aux_phi.matrix().rotated()[nkeep]=qphi.matrix().rotated()[i];
      	      nkeep = nkeep + 1;
      	    }
      	  }
      	  nbase = nkeep ;
      	}else{
      	  nbase = nkeep ;
      	}
	
      }	
    }
  }
}

#ifdef INQ_EIGENSOLVERS_DAVIDSON_UNIT_TEST
#undef INQ_EIGENSOLVERS_DAVIDSON_UNIT_TEST

#include <basis/trivial.hpp>
#include <ions/unitcell.hpp>
#include <operations/matrix_operator.hpp>

#include <catch2/catch.hpp>

TEST_CASE("eigensolvers::davidson", "[eigensolvers::davidson]") {

  using namespace inq;
  using namespace Catch::literals;
	
  const int npoint = 100;
  const int nvec = 4;
  
  basis::trivial bas(npoint);
	
  math::array<complex, 2> identity_matrix({npoint, npoint});
  
  for(int ip = 0; ip < npoint; ip++){
    for(int jp = 0; jp < npoint; jp++){
      identity_matrix[ip][jp] = 0.0;
      if(ip == jp) identity_matrix[ip][jp] = 1.0;
    }
  }
  
  operations::matrix_operator<complex> identity(std::move(identity_matrix));

  SECTION("Diagonal matrix complex"){
  
    math::array<complex, 2> diagonal_matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        diagonal_matrix[ip][jp] = 0.0;
        if(ip == jp) diagonal_matrix[ip][jp] = ip + 1.0;
      }
    }
    
    operations::matrix_operator<complex> diagonal_op(std::move(diagonal_matrix));
    
    basis::field_set<basis::trivial, complex> phi(bas, nvec);

    for(int ip = 0; ip < npoint; ip++){
      for(int ivec = 0; ivec < nvec; ivec++){
        phi.matrix()[ip][ivec] = exp(complex(0.0, (ip*ivec)*0.1));
	//phi.matrix()[ip][ivec] = complex(exp(-abs(ivec-ip)*4.0),0.0);
      }
    }

    operations::orthogonalize(phi);

    const int num_iter = 100;
		
    for(int iter = 0; iter < num_iter; iter++){

      eigensolvers::davidson(diagonal_op, identity, phi);
			
      auto residual = diagonal_op(phi);
      auto eigenvalues = operations::overlap_diagonal(phi, residual);
      operations::shift(-1.0, eigenvalues, phi, residual);
      auto normres = operations::overlap_diagonal(residual);
			
      /*
	tfm::format(std::cout, "  Iteration %4d:\n", iter);
				
	for(int ivec = 0; ivec < phi.set_size(); ivec++){
	tfm::format(std::cout, "    state %4d  evalue = %18.12f  res = %15.10e\n", ivec + 1, real(eigenvalues[ivec]), real(normres[ivec]));
	}
      */

      if(num_iter - 1 == iter){

	CHECK(fabs(eigenvalues[0]) == 1.000000001634_a);
	CHECK(fabs(eigenvalues[1]) == 2.000000003689_a);
	CHECK(fabs(eigenvalues[2]) == 3.000000001955_a);
	CHECK(fabs(eigenvalues[3]) == 4.000000001760_a);
//	CHECK(fabs(eigenvalues[4]) == 5.000000002561_a);
//	CHECK(fabs(eigenvalues[5]) == 6.000000003127_a);
//	CHECK(fabs(eigenvalues[6]) == 7.000000002312_a);
//	CHECK(fabs(eigenvalues[7]) == 8.000000000292_a);
//	CHECK(fabs(eigenvalues[8]) == 8.999999999033_a);
//	CHECK(fabs(eigenvalues[9]) == 9.999999998497_a);
//	CHECK(fabs(eigenvalues[10]) == 10.999999998768_a);
//	CHECK(fabs(eigenvalues[11]) == 11.999999998422_a);

      }
			
    }
 	
  }

#if 0
  SECTION("Periodic Laplacian matrix complex"){
		
    math::array<complex, 2> laplacian_matrix({npoint, npoint});
    
    for(int ip = 0; ip < npoint; ip++){
      for(int jp = 0; jp < npoint; jp++){
        laplacian_matrix[ip][jp] = 0.0;
	identity_matrix[ip][jp] = 0.0;
        if(ip == jp) laplacian_matrix[ip][jp] = -1.0;
        if(ip == jp + 1 or ip == jp - 1) laplacian_matrix[ip][jp] = 2.0;
      }
    }
		
    operations::matrix_operator<complex> laplacian(std::move(laplacian_matrix));
    
    basis::field_set<basis::trivial, complex> phi(bas, nvec);

    phi = 0.0;
    for(int ivec = 0; ivec < nvec; ivec++) phi.matrix()[ivec][ivec] = 1.0;

    for(int iter = 0; iter < 100; iter++){
			
      eigensolvers::davidson(laplacian, identity, phi);
			
      auto residual = laplacian(phi);
      auto eigenvalues = operations::overlap_diagonal(phi, residual);
      operations::shift(eigenvalues, phi, residual, -1.0);
      auto normres = operations::overlap_diagonal(residual);
			
      for(int ivec = 0; ivec < phi.set_size(); ivec++){
	tfm::format(std::cout, " state %4d  evalue = %18.12f  res = %5.0e\n", ivec + 1, real(eigenvalues[ivec]), real(normres[ivec]));
      }
    }
  }
#endif 	

}


#endif


#endif
