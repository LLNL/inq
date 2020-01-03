/* -*- indent-tabs-mode: t -*- */

#ifndef EIGENSOLVER_CONJUGATE_GRADIENT
#define EIGENSOLVER_CONJUGATE_GRADIENT

/*
 Copyright (C) 2019 Xavier Andrade

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
#include <math/d3vector.hpp>
#include <math/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/shift.hpp>

namespace eigensolver {

	template <class operator_type, class preconditioner_type, class field_set_type>
	void conjugate_gradient(const operator_type & ham, const preconditioner_type & prec, basis::field_set<basis::fourier_space, field_set_type> & phi_all){
    
		const int num_steps = 25;
    
    for(int ist = 0; ist < phi_all.set_size(); ist++){
      
      basis::field_set<basis::fourier_space, field_set_type> phi(phi_all.basis(), 1);
      
      phi.matrix().rotated()[0] = phi_all.matrix().rotated()[ist];

      auto hphi = ham(phi);

      auto eigenvalue = operations::overlap_diagonal(hphi, phi);
      auto old_energy = eigenvalue;

      basis::field_set<basis::fourier_space, field_set_type> cg(phi_all.basis(), 1);
      
      for(int istep = 0; istep < num_steps; istep++){

        basis::field_set<basis::fourier_space, field_set_type> g(phi_all.basis(), 1);

        for(long ip = 0; ip < g.basis().size(); ip++) g.matrix()[ip][0] = hphi.matrix()[ip][0] - eigenvalue[0]*phi.matrix()[ip][0];

        auto g0 = g;
        
        prec(g0);

        auto dot = operations::overlap_diagonal(phi, g0);

        //TODO: orthogonalize gg against other states
        
        for(long ip = 0; ip < g.basis().size(); ip++) g0.matrix()[ip][0] = g0.matrix()[ip][0] - dot[0]*phi.matrix()[ip][0];

        auto gg = operations::overlap_diagonal(g0, g);

        std::cout << istep << '\t' << fabs(gg[0]) << std::endl;

        auto gg0 = gg[0];
        
        if(istep == 0){
          cg = g0;          
        } else {
          auto gamma = gg[0]/gg0;
          for(long ip = 0; ip < cg.basis().size(); ip++) cg.matrix()[ip][0] = gamma*cg.matrix()[ip][0] + g0.matrix()[ip][0];

          auto norma = operations::overlap_diagonal(phi, cg);
          for(long ip = 0; ip < cg.basis().size(); ip++) cg.matrix()[ip][0] = cg.matrix()[ip][0] - norma[0]*phi.matrix()[ip][0];
        }

        //cg now contains the conjugate gradient

        auto hcg = ham(cg);

        auto a0 = operations::overlap_diagonal(phi, hcg)[0];
        auto b0 = operations::overlap_diagonal(cg, hcg)[0];
        auto cg0 = operations::overlap_diagonal(cg, cg)[0];
        
        cg0 = sqrt(cg0);

        a0 = 2.0*a0/cg0;
        b0 = b0/(cg0*cg0);
        auto alpha = 2.0*real(eigenvalue[0] - b0);
        auto beta = real(a0)*2.0;        

        auto theta = atan(beta/alpha)*0.5;
        auto stheta = sin(theta);
        auto ctheta = cos(theta);
        auto es1 = alpha*(0.5 - stheta*stheta) + beta*2.0*stheta*ctheta;
        auto stheta2 = sin(theta + M_PI*0.5);
        auto ctheta2 = cos(theta + M_PI*0.5);
        auto es2 = alpha*(0.5 - stheta2*stheta2) + beta*2.0*stheta2*ctheta2;

        if(real(es2) < real(es1)){
          theta += M_PI*0.5;
          a0 = ctheta2;
          b0 = stheta2/cg0;
        } else {
          a0 = ctheta;
          b0 = stheta/cg0;
        }

        for(long ip = 0; ip < cg.basis().size(); ip++){
          phi.matrix()[ip][0] = a0*phi.matrix()[ip][0] - b0*cg.matrix()[ip][0];
          hphi.matrix()[ip][0] = a0*hphi.matrix()[ip][0] - b0*hcg.matrix()[ip][0];
        }
        
      }

      phi_all.matrix().rotated()[ist] = phi.matrix().rotated()[0];
      
    }

  }
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("eigensolver::conjugate_gradient", "[conjugate_gradient]") {


}


#endif


#endif
