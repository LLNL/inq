/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__EIGENSOLVERS__CONJUGATE_GRADIENT
#define INQ__EIGENSOLVERS__CONJUGATE_GRADIENT

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
#include <math/vector3.hpp>
#include <math/array.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <operations/orthogonalize.hpp>
#include <operations/overlap_diagonal.hpp>

#include <caliper/cali.h>

namespace inq {
namespace eigensolvers {

template <class operator_type, class preconditioner_type, class field_set_type>
void conjugate_gradient(const operator_type & ham, const preconditioner_type & prec, basis::field_set<basis::fourier_space, field_set_type> & phi_all){

	CALI_CXX_MARK_FUNCTION;
		
	const int num_iter = 30;
	const double tol = 1.0e-8;
	const double energy_change_threshold = 0.1;
    
	for(int ist = 0; ist < phi_all.set_size(); ist++){

		basis::field_set<basis::fourier_space, field_set_type> phi(phi_all.basis(), 1, phi_all.full_comm());
      
		phi.matrix().rotated()[0] = phi_all.matrix().rotated()[ist];

		operations::orthogonalize_single(phi, phi_all, ist);
      
		auto hphi = ham(phi);

		auto eigenvalue = operations::overlap_diagonal(hphi, phi)[0];
		auto old_energy = eigenvalue;

		double first_delta_e = 0.0;
			
		basis::field_set<basis::fourier_space, field_set_type> cg(phi_all.basis(), 1, phi_all.full_comm());

		complex gg0 = 1.0;
      
		for(int iter = 0; iter < num_iter + 1; iter++){

			eigenvalue = operations::overlap_diagonal(phi, hphi)[0];
        
			basis::field_set<basis::fourier_space, field_set_type> g(phi_all.basis(), 1, phi_all.full_comm());

			auto gm = begin(g.matrix());
			auto phim = begin(phi.matrix());
			auto hphim = begin(hphi.matrix());
				
			gpu::run(g.basis().local_size(),
							 [=] GPU_LAMBDA (auto ip){
								 gm[ip][0] = hphim[ip][0] - eigenvalue*phim[ip][0];
							 });
			
			double res = fabs(operations::overlap_diagonal(g)[0]);

			auto g0 = g;
        
			prec(g0);

			operations::orthogonalize_single(g0, phi_all, ist);

			auto dot = operations::overlap_diagonal(phi, g0);

			gpu::run(g.basis().local_size(),
							 [dd = dot[0],
								g0m = begin(g0.matrix()),
								phim] GPU_LAMBDA (auto ip){
								 g0m[ip][0] = g0m[ip][0] - dd*phim[ip][0];
							 });

			auto gg = operations::overlap_diagonal(g0, g);

			//        std::cout << iter << '\t' << fabs(gg[0]) << std::endl;

			if(sqrt(fabs(gg[0])) < std::numeric_limits<decltype(first_delta_e)>::epsilon()) {
				//					std::cout << "zero gg0" << std::endl;
				break;
			}			
		
			if(iter == 0){
				cg = g0;
				gg0 = gg[0];
			} else {
				auto gamma = gg[0]/gg0;
				gg0 = gg[0];

				auto cgm = begin(cg.matrix());
				auto g0m = begin(g0.matrix());
				
				gpu::run(cg.basis().local_size(),
								 [=] GPU_LAMBDA (auto ip){
									 cgm[ip][0] = gamma*cgm[ip][0] + g0m[ip][0];
								 });
				
				auto norma = operations::overlap_diagonal(phi, cg)[0];

				auto phim = begin(phi.matrix());
				
				gpu::run(cg.basis().local_size(),
								 [=] GPU_LAMBDA (auto ip){
									 cgm[ip][0] = cgm[ip][0] - norma*phim[ip][0];
								 });
									 
			}

			//cg now contains the conjugate gradient

			auto hcg = ham(cg);

			auto a0 = operations::overlap_diagonal(phi, hcg)[0];
			auto b0 = operations::overlap_diagonal(cg, hcg)[0];
			auto cg0 = operations::overlap_diagonal(cg, cg)[0];
        
			cg0 = sqrt(cg0);

			a0 = 2.0*a0/cg0;
			b0 = b0/(cg0*cg0);
			auto alpha = 2.0*real(eigenvalue - b0);
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

			gpu::run(cg.basis().local_size(),
							 [a0, b0,
								phim = begin(phi.matrix()),
								cgm = begin(cg.matrix()),
								hphim = begin(hphi.matrix()),
								hcgm = begin(hcg.matrix())] GPU_LAMBDA (auto ip){
								 phim[ip][0] = a0*phim[ip][0] + b0*cgm[ip][0];
								 hphim[ip][0] = a0*hphim[ip][0] + b0*hcgm[ip][0];
							 });

			//calculate the eigenvalue, this is duplicated
			eigenvalue = operations::overlap_diagonal(phi, hphi)[0];
        
			basis::field_set<basis::fourier_space, field_set_type> g2(phi_all.basis(), 1, phi_all.full_comm());

			gpu::run(g.basis().local_size(),
							 [eigenvalue,
								g2m = begin(g2.matrix()),
								hphim = begin(hphi.matrix()),
								phim = begin(phi.matrix())] GPU_LAMBDA (auto ip){
								 g2m[ip][0] = hphim[ip][0] - eigenvalue*phim[ip][0];
							 });

			res = fabs(operations::overlap_diagonal(g2)[0]);
				
			if(iter > 0){
				if(fabs(eigenvalue - old_energy) < first_delta_e*energy_change_threshold) {
					//						std::cout << "energy_change_threshold " << iter << std::endl;
					break;
				}
			}	else {
				first_delta_e = fabs(eigenvalue - old_energy);
				if(first_delta_e <= 2.0*std::numeric_limits<decltype(first_delta_e)>::epsilon()) {
					//						std::cout << "zero first_delta_e" << std::endl;
					break;
				}
			}

			if(res < tol or iter == num_iter){
				//					std::cout << "res < tol " << ist << '\t' << iter << '\t' << real(eigenvalue[0]) << '\t' << res << std::endl;
				break;
			}
			
			old_energy = eigenvalue;
				
		} // end iteration

		operations::orthogonalize_single(phi, phi_all, ist);

		//normalize
		auto nrm = sqrt(operations::overlap_diagonal(phi)[0]);

		gpu::run(cg.basis().local_size(),
						 [nrm, phim = begin(phi.matrix())] GPU_LAMBDA (auto ip) {
							 phim[ip][0] = phim[ip][0]/nrm;
						 });
		
		// save the newly calculated state
		phi_all.matrix().rotated()[ist] = phi.matrix().rotated()[0];
      
	} // end loop over states

}

}
}

#ifdef INQ_EIGENSOLVERS_CONJUGATE_GRADIENT_UNIT_TEST
#undef INQ_EIGENSOLVERS_CONJUGATE_GRADIENT_UNIT_TEST

#include <catch2/catch.hpp>

#endif

#endif
