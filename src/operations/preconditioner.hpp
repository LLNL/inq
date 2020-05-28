/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__PRECONDITIONER
#define INQ__OPERATIONS__PRECONDITIONER

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

namespace inq {
namespace operations {

class preconditioner {
	
	// Implements the preconditioner of Teter, Payne and Allan, Phys. Rev. B, 40 12255 (1989)
	
public:
	
	GPU_FUNCTION static auto k_function(double x) {

		// The original function of the TPA paper
		auto num = ((8.0*x + 12.0)*x + 18.0)*x + 27.0;
		auto den = (((16.0*x + 8.0)*x + 12.0)*x + 18.0)*x + 27.0;

		/* 
		// A better version by Zhou, Chelikowsky, Gao and Zhou, Commun. Comput. Phys. 18 167 (2015)

		// I won't use it for the moment, since it doesn't seem to
		// change the results and it might have numerical issues. XA

		auto num = ((((32.0*x + 48.0)*x + 72.0)*x + 108.0)*x + 162)*x + 243.0;
		auto den = (((((64.0*x + 32.0)*x + 48.0)*x + 72.0)*x + 108.0)*x + 162)*x + 243.0;
			
		*/
			
		return num/den;
	}
		
	template <class type>
	void operator()(basis::field_set<basis::fourier_space, type> & phi) const {

		math::array<double, 1> expect(phi.set_size(), 0.0);
		math::array<double, 1> norm(phi.set_size(), 0.0);
			
		//calculate the expectation value of the kinetic energy
		//DATAOPERATIONS LOOP + GPU::RUN 4D REDUCTIONS
#ifdef HAVE_CUDA
		gpu::run(phi.set_size(),
						 [expc = begin(expect), nrm = begin(norm), phcub = begin(phi.cubic()), bas = phi.basis()]
						 __device__ (auto ist){
							 for(int ix = 0; ix < bas.sizes()[0]; ix++){
								 for(int iy = 0; iy < bas.sizes()[1]; iy++){
									 for(int iz = 0; iz < bas.sizes()[2]; iz++){
										 auto lapl = -0.5*(-bas.g2(ix, iy, iz));
										 auto phiphi = fabs(phcub[ix][iy][iz][ist]);
										 expc[ist] += lapl*phiphi;
										 nrm[ist] += phiphi;
									 }
								 }
							 }
						 });
#else
		for(int ix = 0; ix < phi.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < phi.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < phi.basis().sizes()[2]; iz++){
					auto lapl = -0.5*(-phi.basis().g2(ix, iy, iz));
					for(int ist = 0; ist < phi.set_size(); ist++){
						auto phiphi = fabs(phi.cubic()[ix][iy][iz][ist]);
						expect[ist] += lapl*phiphi;
						norm[ist] += phiphi;
					}
				}
			}
		}
#endif

		//REDUCE GRID expect norm

		//DATAOPERATIONS LOOP + GPU::RUN 4D
#ifdef HAVE_CUDA
		gpu::run(phi.set_size(), phi.basis().sizes()[2], phi.basis().sizes()[1], phi.basis().sizes()[0], 
						 [expc = begin(expect), nrm = begin(norm), phcub = begin(phi.cubic()), bas = phi.basis()] __device__
						 (auto ist, auto iz, auto iy, auto ix){
							 auto lapl = -0.5*(-bas.g2(ix, iy, iz));
							 phcub[ix][iy][iz][ist] = k_function(lapl*nrm[ist]/expc[ist])*phcub[ix][iy][iz][ist];
						 });
#else
		for(int ix = 0; ix < phi.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < phi.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < phi.basis().sizes()[2]; iz++){
						
					auto lapl = -0.5*(-phi.basis().g2(ix, iy, iz));
					for(int ist = 0; ist < phi.set_size(); ist++) phi.cubic()[ix][iy][iz][ist] *= k_function(lapl*norm[ist]/expect[ist]);
						
				}
			}
		}
#endif

	}

	template <class type>
	void operator()(basis::field_set<basis::real_space, type> & phi) const {
			
		auto fphi = operations::space::to_fourier(phi);
		operator()(fphi);
		phi = operations::space::to_real(fphi);

	}


private:

};
	
}
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::precondition", "[precondition]") {

	using namespace inq;
	using namespace Catch::literals;

}


#endif

#endif
